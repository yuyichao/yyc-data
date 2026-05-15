ion1 = -1
ion2 = 1
nions = 13
nseg = 50

# Optimizer settings
pitime = 3.5  # us, corresponds to max rabi frequency per base function
τmin = 50 / nseg    # min segment length in us
τmax = 150 / nseg   # max segment length in us
maxtime = 10  # max seconds to run optimization
min_mode_index = 1  # lower bound for detune during gate
max_mode_index = min_mode_index + 1  # max bound for detune during gate
params_file = "072125_goldparams_13ions.json"
println("Gate time between $(τmin * nseg) and $(τmax * nseg) μs")
println("Lower bound on pi time required: $(pitime) μs")

using GoldGates
using MSSim: Optimizers as Opts, SegSeq as SS, SymLinear as SL, Sequence as Seq, Utils as U
using NLopt
using Statistics
using JSON
using StaticArrays

function get_am_cbs(NSeg)
    return ntuple(i -> begin
                      prev = (i - 1) / (NSeg / 2) - 1
                      mid = i / (NSeg / 2) - 1
                      next = (i + 1) / (NSeg / 2) - 1
                      function (x)
                          if x < prev || x > next
                              return 0.0
                          elseif x < mid
                              return (x - prev) / (mid - prev)
                          else
                              return (next - x) / (next - mid)
                          end
                      end
                  end, NSeg - 1)
end

# Objective function for optimization
function _objfunc(vals)
    NModes = length(vals) ÷ 3
    @assert length(vals) == NModes * 3 + 1

    dis = sum(vals[1:NModes])
    disδ = sum(vals[NModes + 1:NModes * 2])
    area = sum(abs.(vals[NModes * 2 + 1:NModes * 3]))
    τ = vals[end]

    return (5 * dis + disδ / 100 + 1e-5) / abs(area) * τ
end

function setup_modes(sysparams, ion1, ion2, nions)
    modes = Seq.Modes()
    idx1 = ion1 + (nions + 1) ÷ 2
    idx2 = ion2 + (nions + 1) ÷ 2
    for i in 1:nions
        push!(modes,
            2π * sysparams.modes.radial1[i],
            sysparams.participation_factors[i][idx1] * sysparams.participation_factors[i][idx2] * sysparams.lamb_dicke_parameters[i]^2)
    end
    return modes
end

struct AvgAreaObjCallback{NModes}
    dis_weight::Float64
    disδ_weight::Float64
    area_weights::MVector{NModes,Float64}
    function AvgAreaObjCallback(NModes, dis_weight, disδ_weight, area_weights)
        cb = new{NModes}(dis_weight, disδ_weight, MVector{NModes,Float64}(undef))
        cb.area_weights .= area_weights
        return cb
    end
end

function (obj::AvgAreaObjCallback{NModes})(vals) where NModes
    @assert length(vals) == NModes * 3 + 1
    weights = obj.area_weights
    @inbounds begin
        v1 = muladd(vals[1], obj.dis_weight, vals[1 + NModes] * obj.disδ_weight)
        v2 = abs(vals[1 + NModes * 2]) * weights[1]
    end
    @inbounds @simd ivdep for i in 2:NModes
        v1 = muladd(vals[i], obj.dis_weight, muladd(vals[i + NModes], obj.disδ_weight, v1))
        v2 = muladd(abs(vals[i + NModes * 2]), weights[i], v2)
    end
    v1 = v1 + 1e-5
    return v1 / v2 * vals[end]
end

function avg_area_obj(nseg, modes, pmask;
                      freq=Seq.FreqSpec(), amp=Seq.AmpSpec(),
                      dis_weight=5, disδ_weight=0.01)
    nmodes = length(modes.modes)
    dis_args = ntuple(i->(:dis2, i), nmodes)
    disδ_args = ntuple(i->(:disδ2, i), nmodes)
    area_args = ntuple(i->(:area, i), nmodes)

    mask_dis_area = SS.ValueMask(true, true, true, false, true, false)
    buf = SL.ComputeBuffer{nseg,Float64}(Val(SS.mask_allδ), Val(SS.mask_allδ))

    area_weights = zeros(nmodes)
    area_weights[1] = 1
    area_weights[2] = 1
    return Seq.Objective(pmask, (dis_args..., disδ_args..., area_args..., (:τ, 0)),
                         Opts.autodiff(AvgAreaObjCallback(nmodes, dis_weight,
                                                          disδ_weight, area_weights)),
                         modes, buf, freq=freq, amp=amp)
end

function setup_model(nseg, modes)
    return avg_area_obj(nseg, modes, SL.pmask_full,
                        freq=Seq.FreqSpec(false, sym=true),
                        amp=Seq.AmpSpec(cb=get_am_cbs(nseg), sym=false))
end

function setup_optimizer(nlmodel, sysparams; pitime=30, τmin=5, τmax=50, maxtime=10, min_mode_index=1, max_mode_index=3)
    Ωmax = π / (2 * pitime)
    ωmin = 2π * sysparams.modes.radial1[min_mode_index]
    ωmax = 2π * sysparams.modes.radial1[max_mode_index]

    nargs = Seq.nparams(nlmodel)
    tracker = Opts.NLVarTracker(nargs)
    for Ω in nlmodel.param.Ωs
        Opts.set_bound!(tracker, Ω, -Ωmax, Ωmax)
    end
    Opts.set_bound!(tracker, nlmodel.param.τ, τmin, τmax)
    for ω in nlmodel.param.ωs
        Opts.set_bound!(tracker, ω, ωmin, ωmax)
    end

    opt = NLopt.Opt(:LD_LBFGS, nargs)
    NLopt.min_objective!(opt, nlmodel)
    NLopt.lower_bounds!(opt, Opts.lower_bounds(tracker))
    NLopt.upper_bounds!(opt, Opts.upper_bounds(tracker))
    NLopt.maxtime!(opt, maxtime)

    return opt, tracker
end

function perturb_params(params, tracker, scale)
    result = copy(params)
    lb = Opts.lower_bounds(tracker)
    ub = Opts.upper_bounds(tracker)
    for i in eachindex(result)
        r = (ub[i] - lb[i]) * scale
        result[i] = clamp(params[i] + randn() * r, lb[i], ub[i])
    end
    return result
end

function run_optimization!(opt, tracker, nlmodel; threshold=-Inf,
                           initial_params=nothing, perturbation=0.05)
    iterations = initial_params === nothing ? 500 : 50
    best_obj = Inf
    best_params = nothing

    @time for i in 1:iterations
        if initial_params === nothing
            start_params = Opts.init_vars!(tracker)
        elseif i == 1
            start_params = copy(initial_params)
        else
            start_params = perturb_params(initial_params, tracker, perturbation)
        end
        obj, params, ret = NLopt.optimize(opt, start_params)
        if getfield(NLopt, ret) < 0
            continue
        end
        if obj < best_obj
            best_obj = obj
            area = nlmodel(Val((:area, 0)), params)
            best_status = (
                obj = obj,
                dis = nlmodel(Val((:dis2, 0)), params),
                disδ = nlmodel(Val((:disδ2, 0)), params),
                area = area,
                areaδ = nlmodel(Val((:areaδ, 0)), params),
                total_t = nlmodel(Val((:τ, 0)), params),
                # Ωmax = sum(params[Ω] for Ω in nlmodel.param.Ωs),
                Ωmax = maximum(abs(params[Ω]) for Ω in nlmodel.param.Ωs),
            )
            println(best_status)
            best_params = params
        end
        if best_obj < threshold
            break
        end
    end

    return best_params, best_obj
end

sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
end

modes = setup_modes(sysparams, ion1, ion2, nions)
nlmodel = setup_model(nseg, modes)
opt, tracker = setup_optimizer(nlmodel, sysparams; pitime, τmin, τmax, maxtime, min_mode_index, max_mode_index)

best_params, best_obj = run_optimization!(opt, tracker, nlmodel)

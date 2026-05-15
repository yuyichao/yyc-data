ion1 = -1
ion2 = 1
nions = 13
nseg = 30

# Optimizer settings
pitime = 10  # us, corresponds to max rabi frequency per base function
τmin = 3    # min segment length in us
τmax = 10   # max segment length in us
maxtime = 10  # max seconds to run optimization
min_mode_index = 1  # lower bound for detune during gate
max_mode_index = min_mode_index + 1  # max bound for detune during gate
params_file = "072125_goldparams_13ions.json"
println("Gate time between $(τmin * nseg) and $(τmax * nseg) μs")
println("Lower bound on pi time required: $(pitime/cld(nseg+1, 2)) μs")

using GoldGates
using MSSim: Optimizers as Opts, SegSeq as SS, SymLinear as SL, Sequence as Seq, Utils as U
using NLopt
using Statistics
using JSON

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

function amp_base_funcs(n::Integer; atol::Real = 1e-12)
    @assert n ≥ 1 "n must be at least 1"
    m = n + 1
    grid = [((i - 1) / (n / 2) - 1) for i in 1:m]
    fs_vec = [
        let k = k, grid = grid, m = m, atol = atol
            x -> begin
                exclude = k - 1
                lo = 1 + exclude
                hi = m - exclude
                if lo > hi
                    return 0.0
                end
                for j in lo:hi
                    if isapprox(x, grid[j]; atol = atol)
                        return 1.0
                    end
                end
                return 0.0
            end
        end for k in 1:cld(n+1, 2)
    ]
    return tuple(fs_vec...)
end

# Objective function for optimization
function _objfunc(vals)
    dis = vals[1]
    disδ = vals[2]
    area = vals[3]
    areaδ = vals[4]
    τ = vals[5]

    return 5 * dis + disδ / 100 + (abs(area) - π / 2)^2 * 100
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

function setup_model(nseg, modes)
    objfunc = Opts.autodiff(_objfunc)
    buf_opt = SL.ComputeBuffer{nseg,Float64}(Val(SS.mask_allδ), Val(SS.mask_allδ))
    # amp_funcs = amp_base_funcs(nseg)
    amp_funcs = get_am_cbs(nseg)
    nlmodel = Seq.Objective(SL.pmask_full,
        ((:dis2, 0), (:disδ2, 0), (:area, 0),
            (:areaδ, 0), (:τ, 0)),
        objfunc, modes, buf_opt,
        freq=Seq.FreqSpec(false, sym=true),
        amp=Seq.AmpSpec(cb=amp_funcs, sym=true))
    return nlmodel
end

function setup_optimizer(nlmodel, sysparams; pitime=30, τmin=5, τmax=50, maxtime=10, min_mode_index=1, max_mode_index=3)
    Ωmax = π / (2 * pitime) * 3
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
                areaε = abs(area) - π / 2,
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

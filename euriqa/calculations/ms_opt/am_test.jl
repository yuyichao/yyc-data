#

include("am_shared.jl")

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
using NLopt
using Statistics
using JSON

function setup_modes(sysparams, ion1, ion2, nions)
    modes = Seq.Modes()
    idx1 = ion1 + (nions + 1) ÷ 2
    idx2 = ion2 + (nions + 1) ÷ 2
    for i in 1:nions
        push!(modes, 2π * sysparams.modes.radial1[i])
    end
    return modes
end

function setup_model(nseg, modes)
    return avg_area_obj(nseg, modes, SL.pmask_full,
                        freq=Seq.FreqSpec(false, sym=true),
                        amp=Seq.AmpSpec(cb=get_am_cbs(nseg), sym=false))
end

function setup_optimizer(nlmodel, sysparams; τmin=5, τmax=50, maxtime=10, min_mode_index=1, max_mode_index=3)
    ωmin = 2π * sysparams.modes.radial1[min_mode_index]
    ωmax = 2π * sysparams.modes.radial1[max_mode_index]

    nargs = Seq.nparams(nlmodel)
    tracker = Opts.NLVarTracker(nargs)
    for Ω in nlmodel.param.Ωs
        Opts.set_bound!(tracker, Ω, -1, 1)
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

function run_optimization!(opt, tracker, nlmodel; threshold=-Inf,
                           initial_params=nothing, perturbation=0.05)
    iterations = 100
    best_obj = Inf
    best_params = nothing

    @time for i in 1:iterations
        start_params = Opts.init_vars!(tracker)
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
opt, tracker = setup_optimizer(nlmodel, sysparams; τmin, τmax, maxtime, min_mode_index, max_mode_index)

best_params, best_obj = run_optimization!(opt, tracker, nlmodel)

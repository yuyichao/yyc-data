#

include("am_shared.jl")


using NLopt

using AMO.Utils: ThreadObjectPool, eachobj

using Base.Threads

using GoldGates

struct Candidate
    param::Vector{Float64}
    props::Union{Nothing,Seq.SolutionProperties}
end

struct PreOptimizer{NSeg,PreObj,Sum}
    ωs::Vector{Float64}
    modes::Seq.Modes

    tmin::Float64
    tmax::Float64
    ntimes::Int

    pre_obj::PreObj
    pre_tracker::Opts.NLVarTracker
    pre_opt::NLopt.Opt

    sum::Sum
    args_buff::Vector{Float64}
    rawparams_buff::Vector{Float64}
    candidates::Vector{Candidate}

    function PreOptimizer{NSeg}(ωs; tmin, tmax, ntimes=11, ωmin, ωmax,
                                maxiter=2500, disδ_weight=0.01) where NSeg
        nions = length(ωs)
        modes = Seq.Modes()
        for ω in ωs
            push!(modes, ω)
        end

        freq_spec = Seq.FreqSpec(false, sym=true)
        amp_spec = Seq.AmpSpec(cb=get_am_cbs(NSeg), sym=false)

        pre_obj = avg_area_obj(NSeg, modes, SL.pmask_full,
                               freq=freq_spec, amp=amp_spec, disδ_weight=disδ_weight)

        nargs = Seq.nparams(pre_obj)
        pre_tracker = Opts.NLVarTracker(nargs)
        for Ω in pre_obj.param.Ωs
            Opts.set_bound!(pre_tracker, Ω, -1, 1)
        end
        @assert length(pre_obj.param.ωs) == 1
        for ω in pre_obj.param.ωs
            Opts.set_bound!(pre_tracker, ω, ωmin, ωmax)
        end

        pre_opt = NLopt.Opt(:LD_LBFGS, nargs)
        precompile(pre_obj, (Vector{Float64}, Vector{Float64}))
        NLopt.min_objective!(pre_opt, pre_obj)
        NLopt.maxeval!(pre_opt, maxiter)

        s = Seq.Summarizer{NSeg}()
        return new{NSeg,typeof(pre_obj),typeof(s)}(
            ωs, modes, tmin, tmax, ntimes,
            pre_obj, pre_tracker, pre_opt,
            s, Vector{Float64}(undef, nargs), Vector{Float64}(undef, NSeg * 5),
            Candidate[]
        )
    end
end

function update_bounds!(o::PreOptimizer)
    NLopt.lower_bounds!(o.pre_opt, Opts.lower_bounds(o.pre_tracker))
    NLopt.upper_bounds!(o.pre_opt, Opts.upper_bounds(o.pre_tracker))
    return
end

function opt_one!(o::PreOptimizer)
    objval, args, ret = NLopt.optimize!(o.pre_opt,
                                        Opts.init_vars!(o.pre_tracker, o.args_buff))
    if getfield(NLopt, ret)::NLopt.Result < 0
        return false
    end
    raw_params = Seq.RawParams(o.pre_obj, args, buff=o.rawparams_buff)
    props = get(o.sum, raw_params, o.modes)
    nions = length(o.ωs)

    dis = 0.0
    disδ = 0.0
    max_area = 0.0
    @inbounds @simd ivdep for i in nions
        dis += abs2(props.dis[i])
        disδ += abs2(props.disδ[i])
        max_area = max(max_area, abs(props.area[i]))
    end
    @show objval, dis, disδ, max_area
    if dis < 1e-4 * nions && max_area >= 5
        push!(o.candidates, Candidate(copy(args), props))
        return true
    end
    return false
end

function set_time_range!(o::PreOptimizer, τmin, τmax)
    Opts.set_bound!(o.pre_tracker, o.pre_obj.param.τ, τmin, τmax)
    update_bounds!(o)
end

function opt_all_rounds!(@specialize(cb), pool::ThreadObjectPool{PreOpt},
                         nrounds) where {NSeg,PreOpt<:PreOptimizer{NSeg}}
    o0 = get(pool)
    τs = range(o0.tmin / NSeg, o0.tmax / NSeg, o0.ntimes + 1)
    nions = length(o0.ωs)
    ntimes = o0.ntimes
    put!(pool, o0)
    @threads :greedy for time_idx in 1:ntimes
        o = get(pool)
        set_time_range!(o, τs[time_idx], τs[time_idx + 1])
        for _ in 1:nrounds
            opt_one!(o)
        end
        cb(o.candidates)
        put!(pool, o)
    end
end

function opt_all_rounds!(pool::ThreadObjectPool{PreOpt}, nrounds,
                         candidates=Candidate[]) where {PreOpt<:PreOptimizer}
    opt_all_rounds!(_->nothing, pool, nrounds)
    for o in eachobj(pool)
        append!(candidates, o.candidates)
        empty!(o.candidates)
    end
    return candidates
end

params_file = "072125_goldparams_13ions.json"
sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
end

ωs = 2π .* sysparams.modes.radial1

const pre_pool = ThreadObjectPool() do
    return PreOptimizer{50}(ωs;
                            tmin=50, tmax=150, ntimes=1,
                            ωmin=(ωs[1] + ωs[2]) / 2, ωmax=(ωs[1] + ωs[2]) / 2)
end
candidates = @time opt_all_rounds!(pre_pool, 1000, Candidate[])
@show length(candidates)

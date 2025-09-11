#!/usr/bin/julia

import MSSim.Optimizers as Opts
import MSSim.SegSeq as SS
import MSSim.SymLinear as SL
import MSSim.Sequence as Seq
import MSSim.Utils as U

using JSON
using NLopt
using Base.Threads
using Printf

using GoldGates.ThreadUtils: ThreadObjectPool

function set_mode_weight!(weights, ηs, bij, ion1, ion2)
    nions = length(ηs)
    for i in 1:nions
        weights[i] = bij[i, ion1] * bij[i, ion2] * ηs[i]^2
    end
    return weights
end

include("candidate_pb.jl")

Base.Dict(c::Candidate) = Dict("param"=>c.param, "props"=>Dict(c.props))
Candidate(d::Dict{<:AbstractString}) = Candidate(copy(d["param"]),
                                                 Seq.SolutionProperties(d["props"]))

function load_candidates_file(io::IO)
    data = JSON.parse(io)
    meta = get(data, "meta", nothing)
    return meta, Candidate.(data["candidates"])
end

function load_candidates_pb(io::IO)
    decoder = PB.ProtoDecoder(io)
    candidates = PB.decode(decoder, Candidates)
    if candidates.meta == ""
        meta = nothing
    else
        meta = JSON.parse(candidates.meta)
    end
    return meta, candidates.candidates
end

function load_candidates_files(files; candidates=Candidate[], meta=nothing)
    results = Dict{String,Tuple{Dict,Vector{Candidate}}}()
    lock = ReentrantLock()
    @threads :greedy for f in files
        if endswith(f, ".binpb")
            res = open(load_candidates_pb, f)
        else
            res = open(load_candidates_file, f)
        end
        @lock lock results[f] = res
    end
    for f in files
        file_meta, file_candidates = results[f]
        if meta === nothing
            meta = file_meta
        elseif file_meta != meta
            error("Metadata mismatch")
        end
        append!(candidates, file_candidates)
    end
    return meta, candidates
end

function load_candidates_dir(dir; kwargs...)
    return load_candidates_files(readdir(dir, join=true); kwargs...)
end

function load_candidates_dirs(dirs; kwargs...)
    files = String[]
    for dir in dirs
        append!(files, readdir(dir, join=true))
    end
    return load_candidates_files(files; kwargs...)
end

function save_candidates(prefix, candidates, meta; block_size=2000, use_pb=false)
    ncandidates = length(candidates)
    meta_str = nothing
    if use_pb && meta !== nothing
        meta_str = JSON.json(meta)
    end
    @threads :greedy for (i, start_idx) in enumerate(1:block_size:ncandidates)
        end_idx = min(start_idx + block_size - 1, ncandidates)
        if use_pb
            open("$(prefix)$(@sprintf("%06d", i)).binpb", "w") do io
                encoder = PB.ProtoEncoder(io)
                PB.encode(encoder, Candidates(candidates[start_idx:end_idx], meta_str))
            end
        else
            open("$(prefix)$(@sprintf("%06d", i)).json", "w") do io
                d = Dict("meta"=>meta,
                         "candidates"=>Dict.(@view candidates[start_idx:end_idx]))
                JSON.print(io, d, 2)
            end
        end
    end
end

struct PreOptimizer{NSeg,PreObj,Sum}
    ωs::Vector{Float64}
    ωs2::Vector{Float64}
    modes::Seq.Modes
    modes2::Seq.Modes

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

    function PreOptimizer{NSeg}(ωs, ωs2; tmin, tmax, ntimes=11, ωmin, ωmax,
                                amp_ratio=0.7, maxiter=2500) where NSeg
        nions = length(ωs)
        modes = Seq.Modes()
        for ω in ωs
            push!(modes, ω)
        end
        modes2 = Seq.Modes()
        for ω in ωs
            push!(modes2, ω)
        end
        for ω in ωs2
            push!(modes2, ω)
        end

        freq_spec = Seq.FreqSpec(true, sym=false)
        amp_spec = Seq.AmpSpec(cb=U.BlackmanStartEnd{amp_ratio}())

        pre_obj = Opts.abs_area_obj(NSeg, modes, SL.pmask_tfm,
                                    freq=freq_spec, amp=amp_spec,
                                    dis_weights=fill(1, nions),
                                    disδ_weights=fill(0.01, nions),
                                    area_weights=zeros(nions))
        nargs = Seq.nparams(pre_obj)
        pre_tracker = Opts.NLVarTracker(nargs)
        Opts.set_bound!(pre_tracker, pre_obj.param.Ωs[1], 1, 1)
        for ω in pre_obj.param.ωs
            Opts.set_bound!(pre_tracker, ω, ωmin, ωmax)
        end

        pre_opt = NLopt.Opt(:LD_LBFGS, nargs)
        precompile(pre_obj, (Vector{Float64}, Vector{Float64}))
        NLopt.min_objective!(pre_opt, pre_obj)
        NLopt.maxeval!(pre_opt, maxiter)

        s = Seq.Summarizer{NSeg}()
        return new{NSeg,typeof(pre_obj),typeof(s)}(
            ωs, ωs2, modes, modes2, tmin, tmax, ntimes,
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
    props = get(o.sum, raw_params, o.modes2)
    nions = length(o.ωs)

    dis = 0.0
    disδ = 0.0
    max_area = 0.0
    @inbounds @simd ivdep for i in nions
        dis += abs2(props.dis[i])
        disδ += abs2(props.disδ[i])
        max_area = max(max_area, abs(props.area[i]))
    end
    if dis < 1e-6 * nions && disδ < 1e-4 * nions && max_area >= 100
        push!(o.candidates, Candidate(copy(args), props))
        return true
    end
    return false
end

function set_mode!(o::PreOptimizer, mode_idx)
    nions = length(o.ωs)
    @assert 1 <= mode_idx <= nions
    area_weights = o.pre_obj.obj.area_weights
    for i in 1:nions
        area_weights[i] = i == mode_idx
    end
end

function set_time_range!(o::PreOptimizer, τmin, τmax)
    Opts.set_bound!(o.pre_tracker, o.pre_obj.param.τ, τmin, τmax)
    update_bounds!(o)
end

function opt_all_rounds!(pool::ThreadObjectPool{PreOpt},
                         nrounds, candidates) where {NSeg,PreOpt<:PreOptimizer{NSeg}}
    o0 = get(pool)
    τs = range(o0.tmin / NSeg, o0.tmax / NSeg, o0.ntimes + 1)
    nions = length(o0.ωs)
    ntimes = o0.ntimes
    put!(pool, o0)
    @threads :greedy for (mode_idx, time_idx) in Iterators.product(1:nions, 1:ntimes)
        o = get(pool)
        set_mode!(o, mode_idx)
        set_time_range!(o, τs[time_idx], τs[time_idx + 1])
        for _ in 1:nrounds
            opt_one!(o)
        end
        @lock pool.lock begin
            append!(candidates, o.candidates)
        end
        empty!(o.candidates)
        put!(pool, o)
    end
    return candidates
end

struct GateModeInfo
    ωs::Vector{Float64}
    ηs::Vector{Float64}
    bij::Matrix{Float64}
    weights_buff::Vector{Float64}
    GateModeInfo(ωs, ηs, bij) = new(ωs, ηs, bij, Vector{Float64}(undef, length(ωs)))
end

get_weights!(info::GateModeInfo, ion1, ion2) =
    set_mode_weight!(info.weights_buff, info.ηs, info.bij, ion1, ion2)

mutable struct PairChecker
    const weights::Vector{Float64}
    const areaδ::Vector{Float64}
    const minarea::Float64
    candidates_processed::Int
    passed::Bool
    PairChecker(weights, minarea) = new(weights, Float64[], minarea, 0, false)
end

function _check_areaδ(c::PairChecker)
    c1 = 0
    c2 = 0
    for areaδ in c.areaδ
        if areaδ <= 10
            return true
        elseif areaδ <= 20
            c1 += 1
        elseif areaδ <= 40
            c2 += 1
        end
    end
    return c1 >= 10 || (c1 < 6 && c2 >= 30)
end

function check(c::PairChecker, candidates)
    if c.passed
        return true
    end
    nmodes = length(c.weights)
    ncandidates = length(candidates)
    # areas = Float64[]
    # max_area = 0.0
    # max_area_time = 0.0
    for i in c.candidates_processed + 1:ncandidates
        cand = candidates[i]
        area = abs(sum(c.weights[j] * cand.props.area[j] for j in 1:nmodes))
        # if max_area < area
        #     max_area = area
        #     max_area_time = cand.props.total_time
        # end
        # push!(areas, area)
        if area < c.minarea
            continue
        end
        areaδ = abs(sum(c.weights[j] * cand.props.areaδ[j] for j in 1:nmodes)) / area
        push!(c.areaδ, areaδ)
    end
    c.candidates_processed = ncandidates
    if _check_areaδ(c)
        c.passed = true
        return true
    end
    # println("  $(length(c.areaδ)), $(max_area) in $(max_area_time)")
    return false
end

struct PairOptimizer{NSeg,AreaObj,Sum}
    ωs::Vector{Float64}
    ηs::Vector{Float64}
    bij::Matrix{Float64}
    ωs2::Vector{Float64}
    modes::Seq.Modes

    Ωmax::Float64

    area_obj::AreaObj
    area_tracker::Opts.NLVarTracker
    area_opt::NLopt.Opt

    sum::Sum
    weights_buff::Vector{Float64}
    args_buff::Vector{Float64}
    rawparams_buff::Vector{Float64}

    function PairOptimizer{NSeg}(ωs, ηs, bij, ωs2;
                                 tmin, tmax, ntimes=11,
                                 ωmin, ωmax, δω=2π * 0.0005,
                                 Ωmax=0.4, amp_ratio=0.7,
                                 area_maxiter=10000) where NSeg
        nions = length(ωs)
        modes = Seq.Modes()
        for ω in ωs
            push!(modes, ω)
        end
        for ω in ωs
            push!(modes, ω - δω)
        end
        for ω in ωs
            push!(modes, ω + δω)
        end
        for ω in ωs2
            push!(modes, ω)
        end

        freq_spec = Seq.FreqSpec(true, sym=false)
        amp_spec = Seq.AmpSpec(cb=U.BlackmanStartEnd{amp_ratio}())

        dis_weights = [fill(1, nions); fill(0.002, nions * 2); fill(0.0005, nions)]
        disδ_weights = [fill(0.01, nions); fill(0, nions * 2); fill(0.0000001, nions)]
        area_targets = [Opts.AreaTarget(1, area_weights=zeros(nions),
                                        areaδ_weights=zeros(nions)),
                        Opts.AreaTarget(nions + 1, area_weights=zeros(nions)),
                        Opts.AreaTarget(nions * 2 + 1, area_weights=zeros(nions))]

        area_obj = Opts.target_area_obj(NSeg, modes, SL.pmask_full,
                                        freq=freq_spec, amp=amp_spec,
                                        dis_weights=dis_weights,
                                        disδ_weights=disδ_weights,
                                        area_targets=area_targets)
        nargs = Seq.nparams(area_obj)
        area_tracker = Opts.NLVarTracker(nargs)
        Opts.set_bound!(area_tracker, area_obj.param.Ωs[1], 0.01 * Ωmax, Ωmax)
        for ω in area_obj.param.ωs
            Opts.set_bound!(area_tracker, ω, ωmin, ωmax)
        end
        area_opt = NLopt.Opt(:LD_LBFGS, nargs)
        precompile(area_obj, (Vector{Float64}, Vector{Float64}))
        NLopt.min_objective!(area_opt, area_obj)
        NLopt.maxeval!(area_opt, area_maxiter)

        s = Seq.Summarizer{NSeg}()
        return new{NSeg,typeof(area_obj),typeof(s)}(
            ωs, ηs, bij, ωs2, modes, Ωmax,
            area_obj, area_tracker, area_opt,
            s, Vector{Float64}(undef, nions), Vector{Float64}(undef, nargs),
            Vector{Float64}(undef, NSeg * 5)
        )
    end
end

function update_bounds!(o::PairOptimizer)
    NLopt.lower_bounds!(o.area_opt, Opts.lower_bounds(o.area_tracker))
    NLopt.upper_bounds!(o.area_opt, Opts.upper_bounds(o.area_tracker))
    return
end

get_weights!(o::PairOptimizer, ion1, ion2) =
    set_mode_weight!(o.weights_buff, o.ηs, o.bij, ion1, ion2)

function opt_pair!(o::PairOptimizer, args, area, weights)
    args[o.area_obj.param.Ωs[1]] .*= sqrt((π / 2) / area)

    τ0 = args[o.area_obj.param.τ]
    Opts.set_bound!(o.area_tracker, o.area_obj.param.τ, 0.95 * τ0, 1.05 * τ0)
    update_bounds!(o)

    area_tgt1 = o.area_obj.area_targets[1]
    area_tgt2 = o.area_obj.area_targets[2]
    area_tgt3 = o.area_obj.area_targets[3]

    area_tgt1.target = π / 2 * 10
    area_tgt1.area_weights .= weights * 10
    area_tgt1.areaδ_weights .= weights .* 0.01
    area_tgt3.target = π / 2 * 0.001
    area_tgt3.area_weights .= weights .* 0.001
    area_tgt2.target = π / 2 * 0.001
    area_tgt2.area_weights .= weights .* 0.001

    objval, args, ret = NLopt.optimize!(o.area_opt, args)
    if getfield(NLopt, ret)::NLopt.Result < 0
        return
    end
    raw_params = Seq.RawParams(o.area_obj, args, buff=o.rawparams_buff)
    props = get(o.sum, raw_params, o.modes)
    return Candidate(args, props)
end

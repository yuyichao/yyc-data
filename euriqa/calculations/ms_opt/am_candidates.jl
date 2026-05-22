#

include("am_shared.jl")


using NLopt

using AMO.Utils: ThreadObjectPool, eachobj

using Base.Threads

using GoldGates

import ProtoBuf as PB

struct Candidate
    τ::Float64
    ω::Float64
    Ωs::Vector{Float64}
    props::Union{Nothing,Seq.SolutionProperties}
end
function Candidate(args::AbstractVector, param::Seq.ModSpec, props)
    @assert length(param.ωs) == 1
    @assert length(args) == length(param.Ωs) + 2
    return Candidate(args[param.τ], args[param.ωs[1]],
                     [args[Ω] for Ω in param.Ωs], props)
end

struct Candidates
    candidates::Vector{Candidate}
    meta::String
end

module CandidateSerialization

import ..Candidate, ..Candidates

import ProtoBuf as PB
PB.default_values(::Type{Candidate}) = (;τ = 0.0, ω = 0.0, Ωs = Vector{Float64}(),
                                        props = nothing)
PB.field_numbers(::Type{Candidate}) = (;τ=1, ω=2, Ωs=3, props = 4)
Base.:(==)(v1::Candidate, v2::Candidate) = v1.τ == v2.τ && v1.ω == v2.ω && v1.Ωs == v2.Ωs && v1.props == v2.props
Base.hash(v::Candidate, h::UInt) = hash(v.τ, hash(v.ω, hash(v.Ωs, hash(v.props, hash(:Candidate, h)))))

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:Candidate})
    τ = zero(Float64)
    ω = zero(Float64)
    Ωs = PB.BufferedVector{Float64}()
    props = Ref{Union{Nothing,Seq.SolutionProperties}}(nothing)
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            τ = PB.decode(d, Float64)
        elseif field_number == 2
            ω = PB.decode(d, Float64)
        elseif field_number == 3
            PB.decode!(d, wire_type, param)
        elseif field_number == 4
            PB.decode!(d, props)
        else
            PB.skip(d, wire_type)
        end
    end
    return Candidate(τ, ω, Ωs, props[])
end

function PB.encode(e::PB.AbstractProtoEncoder, x::Candidate)
    initpos = position(e.io)
    x.τ !== zero(Float64) && PB.encode(e, 1, x.τ)
    x.ω !== zero(Float64) && PB.encode(e, 2, x.ω)
    !isempty(x.Ωs) && PB.encode(e, 3, x.Ωs)
    !isnothing(x.props) && PB.encode(e, 4, x.props)
    return position(e.io) - initpos
end
function PB._encoded_size(x::Candidate)
    encoded_size = 0
    x.τ !== zero(Float64) && (encoded_size += PB._encoded_size(x.τ, 1))
    x.ω !== zero(Float64) && (encoded_size += PB._encoded_size(x.ω, 2))
    !isempty(x.Ωs) && (encoded_size += PB._encoded_size(x.Ωs, 3))
    !isnothing(x.props) && (encoded_size += PB._encoded_size(x.props, 4))
    return encoded_size
end

PB.default_values(::Type{Candidates}) = (;candidates = Vector{Candidate}(), meta = "")
PB.field_numbers(::Type{Candidates}) = (;candidates = 1, meta = 2)
Base.:(==)(v1::Candidates, v2::Candidates) =
    v1.candidates == v2.candidates && v1.meta == v2.meta
Base.hash(v::Candidates, h::UInt) = hash(v.candidates, hash(v.meta, hash(:Candidates, h)))

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:Candidates})
    candidates = PB.BufferedVector{Candidate}()
    meta = ""
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            PB.decode!(d, candidates)
        elseif field_number == 2
            meta = PB.decode(d, String)
        else
            PB.skip(d, wire_type)
        end
    end
    return Candidates(candidates[], meta)
end

function PB.encode(e::PB.AbstractProtoEncoder, x::Candidates)
    initpos = position(e.io)
    !isempty(x.candidates) && PB.encode(e, 1, x.candidates)
    !isempty(x.meta) && PB.encode(e, 2, x.meta)
    return position(e.io) - initpos
end
function PB._encoded_size(x::Candidates)
    encoded_size = 0
    !isempty(x.candidates) && (encoded_size += PB._encoded_size(x.candidates, 1))
    !isempty(x.meta) && (encoded_size += PB._encoded_size(x.meta, 2))
    return encoded_size
end

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
    if dis < 1e-4 * nions && max_area >= 30
        @show objval, dis, disδ, max_area
        push!(o.candidates, Candidate(args, o.pre_obj.param, props))
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

const prefix = ARGS[1]

params_file = "072125_goldparams_13ions.json"
sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
end

ωs = 2π .* sysparams.modes.radial1

const pre_pool = ThreadObjectPool() do
    return PreOptimizer{50}(ωs;
                            tmin=100, tmax=150, ntimes=1,
                            ωmin=(ωs[1] + ωs[2]) / 2, ωmax=(ωs[1] + ωs[2]) / 2)
end
candidates = @time opt_all_rounds!(pre_pool, 100, Candidate[])
@show length(candidates)

open("$(prefix).binpb", "w") do io
    encoder = PB.ProtoEncoder(io)
    PB.encode(encoder, Candidates(candidates, ""))
end

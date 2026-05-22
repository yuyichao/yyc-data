#

import MSSim: Optimizers as Opts, SegSeq as SS, SymLinear as SL, Sequence as Seq, Utils as U
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

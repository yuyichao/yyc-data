#

import ProtoBuf as PB

struct Candidate
    param::Vector{Float64}
    props::Union{Nothing,Seq.SolutionProperties}
end
PB.default_values(::Type{Candidate}) = (;param = Vector{Float64}(), props = nothing)
PB.field_numbers(::Type{Candidate}) = (;param = 1, props = 2)

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:Candidate})
    param = PB.BufferedVector{Float64}()
    props = Ref{Union{Nothing,Seq.SolutionProperties}}(nothing)
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            PB.decode!(d, wire_type, param)
        elseif field_number == 2
            PB.decode!(d, props)
        else
            PB.skip(d, wire_type)
        end
    end
    return Candidate(param[], props[])
end

function PB.encode(e::PB.AbstractProtoEncoder, x::Candidate)
    initpos = position(e.io)
    !isempty(x.param) && PB.encode(e, 1, x.param)
    !isnothing(x.props) && PB.encode(e, 2, x.props)
    return position(e.io) - initpos
end
function PB._encoded_size(x::Candidate)
    encoded_size = 0
    !isempty(x.param) && (encoded_size += PB._encoded_size(x.param, 1))
    !isnothing(x.props) && (encoded_size += PB._encoded_size(x.props, 2))
    return encoded_size
end

struct Candidates
    candidates::Vector{Candidate}
    meta::String
end
PB.default_values(::Type{Candidates}) = (;candidates = Vector{Candidate}(), meta = "")
PB.field_numbers(::Type{Candidates}) = (;candidates = 1, meta = 2)

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

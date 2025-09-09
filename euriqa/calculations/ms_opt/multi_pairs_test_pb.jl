#

import ProtoBuf as PB

struct SolutionInfo
    params::Vector{Float64}
    total_t::Float64
    dis2::Float64
    disδ2::Float64
    area::Vector{Float64}
    areaδ::Vector{Float64}
end
PB.default_values(::Type{SolutionInfo}) = (;params = Vector{Float64}(), total_t = zero(Float64), dis2 = zero(Float64), disδ2 = zero(Float64), area = Vector{Float64}(), areaδ = Vector{Float64}())
PB.field_numbers(::Type{SolutionInfo}) = (;params = 1, total_t = 2, dis2 = 3, disδ2 = 4, area = 5, areaδ = 6)

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:SolutionInfo})
    params = PB.BufferedVector{Float64}()
    total_t = zero(Float64)
    dis2 = zero(Float64)
    disδ2 = zero(Float64)
    area = PB.BufferedVector{Float64}()
    areaδ = PB.BufferedVector{Float64}()
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            PB.decode!(d, wire_type, params)
        elseif field_number == 2
            total_t = PB.decode(d, Float64)
        elseif field_number == 3
            dis2 = PB.decode(d, Float64)
        elseif field_number == 4
            disδ2 = PB.decode(d, Float64)
        elseif field_number == 5
            PB.decode!(d, wire_type, area)
        elseif field_number == 6
            PB.decode!(d, wire_type, areaδ)
        else
            PB.skip(d, wire_type)
        end
    end
    return SolutionInfo(params[], total_t, dis2, disδ2, area[], areaδ[])
end

function PB.encode(e::PB.AbstractProtoEncoder, x::SolutionInfo)
    initpos = position(e.io)
    !isempty(x.params) && PB.encode(e, 1, x.params)
    x.total_t !== zero(Float64) && PB.encode(e, 2, x.total_t)
    x.dis2 !== zero(Float64) && PB.encode(e, 3, x.dis2)
    x.disδ2 !== zero(Float64) && PB.encode(e, 4, x.disδ2)
    !isempty(x.area) && PB.encode(e, 5, x.area)
    !isempty(x.areaδ) && PB.encode(e, 6, x.areaδ)
    return position(e.io) - initpos
end
function PB._encoded_size(x::SolutionInfo)
    encoded_size = 0
    !isempty(x.params) && (encoded_size += PB._encoded_size(x.params, 1))
    x.total_t !== zero(Float64) && (encoded_size += PB._encoded_size(x.total_t, 2))
    x.dis2 !== zero(Float64) && (encoded_size += PB._encoded_size(x.dis2, 3))
    x.disδ2 !== zero(Float64) && (encoded_size += PB._encoded_size(x.disδ2, 4))
    !isempty(x.area) && (encoded_size += PB._encoded_size(x.area, 5))
    !isempty(x.areaδ) && (encoded_size += PB._encoded_size(x.areaδ, 6))
    return encoded_size
end

struct TimeRangeSolution
    total_t_min::Float64
    total_t_max::Float64
    solutions::Vector{SolutionInfo}
end
PB.default_values(::Type{TimeRangeSolution}) = (;total_t_min = zero(Float64), total_t_max = zero(Float64), solutions = Vector{SolutionInfo}())
PB.field_numbers(::Type{TimeRangeSolution}) = (;total_t_min = 1, total_t_max = 2, solutions = 3)

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:TimeRangeSolution})
    total_t_min = zero(Float64)
    total_t_max = zero(Float64)
    solutions = PB.BufferedVector{SolutionInfo}()
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            total_t_min = PB.decode(d, Float64)
        elseif field_number == 2
            total_t_max = PB.decode(d, Float64)
        elseif field_number == 3
            PB.decode!(d, solutions)
        else
            PB.skip(d, wire_type)
        end
    end
    return TimeRangeSolution(total_t_min, total_t_max, solutions[])
end

function PB.encode(e::PB.AbstractProtoEncoder, x::TimeRangeSolution)
    initpos = position(e.io)
    x.total_t_min !== zero(Float64) && PB.encode(e, 1, x.total_t_min)
    x.total_t_max !== zero(Float64) && PB.encode(e, 2, x.total_t_max)
    !isempty(x.solutions) && PB.encode(e, 3, x.solutions)
    return position(e.io) - initpos
end
function PB._encoded_size(x::TimeRangeSolution)
    encoded_size = 0
    x.total_t_min !== zero(Float64) && (encoded_size += PB._encoded_size(x.total_t_min, 1))
    x.total_t_max !== zero(Float64) && (encoded_size += PB._encoded_size(x.total_t_max, 2))
    !isempty(x.solutions) && (encoded_size += PB._encoded_size(x.solutions, 3))
    return encoded_size
end

struct NSegSolution
    nseg::UInt32
    solutions::Vector{TimeRangeSolution}
end
PB.default_values(::Type{NSegSolution}) = (;nseg = zero(UInt32), solutions = Vector{TimeRangeSolution}())
PB.field_numbers(::Type{NSegSolution}) = (;nseg = 1, solutions = 2)

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:NSegSolution})
    nseg = zero(UInt32)
    solutions = PB.BufferedVector{TimeRangeSolution}()
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            nseg = PB.decode(d, UInt32)
        elseif field_number == 2
            PB.decode!(d, solutions)
        else
            PB.skip(d, wire_type)
        end
    end
    return NSegSolution(nseg, solutions[])
end

function PB.encode(e::PB.AbstractProtoEncoder, x::NSegSolution)
    initpos = position(e.io)
    x.nseg != zero(UInt32) && PB.encode(e, 1, x.nseg)
    !isempty(x.solutions) && PB.encode(e, 2, x.solutions)
    return position(e.io) - initpos
end
function PB._encoded_size(x::NSegSolution)
    encoded_size = 0
    x.nseg != zero(UInt32) && (encoded_size += PB._encoded_size(x.nseg, 1))
    !isempty(x.solutions) && (encoded_size += PB._encoded_size(x.solutions, 2))
    return encoded_size
end

#

include("am_shared.jl")
include("ind-modulation.jl")

using GoldGates
import ProtoBuf as PB

params_file = "072125_goldparams_13ions.json"
sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
end

ωs = 2π .* sysparams.modes.radial1

candidates_file = "data_am/am_candidates_20260522.binpb"
candidates = open(candidates_file) do io
    decoder = PB.ProtoDecoder(io)
    candidates = PB.decode(decoder, Candidates)
    return candidates.candidates
end

function check_candidate(c1, c2)
    if !(0.95 < c1.τ / c2.τ < 1.05)
        return true
    end
    d1 = sum(Ω^2 for Ω in c1.Ωs)
    d2 = sum(Ω^2 for Ω in c2.Ωs)
    d12 = sum(Ω1 * Ω2 for (Ω1, Ω2) in zip(c1.Ωs, c2.Ωs))
    if abs2(d12) > d1 * d2 * 0.8
        return false
    end
    return true
end

function filter_candidates!(candidates)
    i = 1
    while i <= length(candidates)
        newc = candidates[i]
        keep = true
        for j in 1:i - 1
            if !check_candidate(candidates[j], newc)
                keep = false
                break
            end
        end
        if keep
            i = i + 1
        else
            deleteat!(candidates, i)
        end
    end
end

@show length(candidates)
filter_candidates!(candidates)
@show length(candidates)

for c in candidates
    pushfirst!(c.Ωs, 0)
    push!(c.Ωs, 0)
end

function compute_area2_cand(c1, c2, ωm)
    compute_area2_am(c1.τ, c2.τ, c1.ω, ωm, c1.Ωs, c2.Ωs)
end

function area2_matrix(candidates, ωm)
    ncands = length(candidates)
    areas = Matrix{Float64}(undef, ncands, ncands)
    @inbounds for i in 1:ncands
        c1 = candidates[i]
        for j in 1:ncands
            c2 = candidates[j]
            areas[j, i] = compute_area2_am(c1.τ, c2.τ, c1.ω, ωm, c1.Ωs, c2.Ωs)
        end
    end
    @inbounds for i in 1:ncands
        for j in 1:i - 1
            areas[j, i] = areas[i, j] = (areas[j, i] + areas[i, j]) / 2
        end
    end
    return areas
end

for ωm in ωs
    @time area2_matrix(candidates, ωm)
end

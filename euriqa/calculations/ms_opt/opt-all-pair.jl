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

using GoldGates: Candidate,
    load_candidates_files, load_candidates_dir, load_candidates_dirs,
    save_candidates
using GoldGates.Optimizers
using AMO.Utils: ThreadObjectPool

function set_mode_weight!(weights, ηs, bij, ion1, ion2)
    nions = length(ηs)
    for i in 1:nions
        weights[i] = bij[i, ion1] * bij[i, ion2] * ηs[i]^2
    end
    return weights
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

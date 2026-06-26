#

include("am_shared.jl")
include("ind-modulation.jl")

using LinearAlgebra
using GoldGates
import ProtoBuf as PB
using NLopt
using MSSim: Optimizers as Opts

params_file = "072125_goldparams_13ions.json"
sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
end

candidates_file = "data_am2/am_candidates_20260625_4.binpb"
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

function couple_matrix(sysparams, modei, nions)
    return [sysparams.participation_factors[modei][ion1] * sysparams.participation_factors[modei][ion2] * sysparams.lamb_dicke_parameters[modei]^2
            for ion1 in 1:nions, ion2 in 1:nions]
end

struct AreaKernel
    nions::Int
    ncands::Int
    all_areas::Vector{Matrix{Float64}}
    all_couple::Vector{Matrix{Float64}}
    cs′::Vector{Matrix{Float64}}
    function AreaKernel(candidates, sysparams)
        ωs = 2π .* sysparams.modes.radial1
        nions = length(ωs)
        ncands = length(candidates)
        return new(nions, ncands,
                   [area2_matrix(candidates, ωm) for ωm in ωs],
                   [couple_matrix(sysparams, modei, nions) for modei in 1:nions],
                   [Matrix{Float64}(undef, ncands, nions) for modei in 1:nions])
    end
end

function _popolate_cs′(k::AreaKernel, cs)
    @inbounds for modei in 1:k.nions
        mul!(k.cs′[modei], k.all_areas[modei], cs)
    end
end

function pair_area(k::AreaKernel, ion1, ion2, cs)
    res = zero(eltype(cs))
    @inbounds for modei in 1:k.nions
        res = muladd(k.all_couple[modei][ion1, ion2],
                     dot(@view(cs[:, ion1]), @view(k.cs′[modei][:, ion2])), res)
    end
    return res
end

function all_pair_areas(k::AreaKernel, cs)
    nions = k.nions
    _popolate_cs′(k, cs)
    res = zeros(eltype(cs), nions, nions)
    for ion1 in 1:nions
        for ion2 in 1:ion1 - 1
            res[ion1, ion2] = res[ion2, ion1] = pair_area(k, ion1, ion2, cs)
        end
    end
    return res
end

@inbounds function area_target_model(k::AreaKernel, tgt, cs, grad)
    nions = k.nions
    res = zero(eltype(cs))
    _popolate_cs′(k, cs)
    has_grad = !isempty(grad)
    if has_grad
        fill!(grad, zero(eltype(grad)))
    end
    @inbounds for ion1 in 1:nions
        for ion2 in 1:ion1 - 1
            d = pair_area(k, ion1, ion2, cs) - tgt[ion2, ion1]
            res = muladd(d, d, res)
            if has_grad
                for modei in 1:nions
                    s = 2 * d * k.all_couple[modei][ion1, ion2]
                    cs′ = k.cs′[modei]
                    @simd for candi in 1:k.ncands
                        grad[candi, ion1] = muladd(s, cs′[candi, ion2], grad[candi, ion1])
                        grad[candi, ion2] = muladd(s, cs′[candi, ion1], grad[candi, ion2])
                    end
                end
            end
        end
    end
    return res
end

function get_tracker(nargs::Int, Ωmax)
    tracker = Opts.NLVarTracker(nargs)
    for i in 1:nargs
        Opts.set_bound!(tracker, i, -Ωmax, Ωmax)
    end
    return tracker
end
function get_tracker(k::AreaKernel, Ωmax)
    return get_tracker(k.ncands * k.nions, Ωmax)
end

function get_objective(kernel, tgt)
    nions = kernel.nions
    nterms = kernel.ncands
    return function (x, grad)
        res = if isempty(grad)
            area_target_model(kernel, tgt, reshape(x, (nterms, nions)), ())
        else
            area_target_model(kernel, tgt, reshape(x, (nterms, nions)),
                              reshape(grad, (nterms, nions)))
        end
    end
end

function get_opt(k::AreaKernel, tgt, tracker::Opts.NLVarTracker; maxtime=10, xtol=1e-6)
    lb = Opts.lower_bounds(tracker)
    ub = Opts.upper_bounds(tracker)
    opt = NLopt.Opt(:LD_CCSAQ, length(lb))
    NLopt.maxtime!(opt, maxtime)
    NLopt.xtol_rel!(opt, xtol)
    NLopt.lower_bounds!(opt, lb)
    NLopt.upper_bounds!(opt, ub)
    NLopt.min_objective!(opt, get_objective(k, tgt))
    return opt
end

function opt_one!(opt, tracker, args_buff)
    objval, args, ret = NLopt.optimize!(opt, Opts.init_vars!(tracker, args_buff))
    if getfield(NLopt, ret)::NLopt.Result < 0
        return
    end
    return objval, args
end

function opt_rep(opt, tracker, n)
    args_buff = Vector{Float64}(undef, ndims(opt))
    best_val = 1
    best_args = nothing
    for i in 1:n
        res = @time opt_one!(opt, tracker, args_buff)
        if res === nothing
            continue
        end
        val, args = res
        if val < best_val
            @show val
            best_val = val
            best_args = copy(args)
        end
    end
    return best_val, best_args
end

function eval_cand(c, t)
    tr = t / c.τ
    i = floor(Int, tr) + 1
    if i < 1 || i >= length(c.Ωs)
        return 0.0
    end
    Ω1 = c.Ωs[i]
    Ω2 = c.Ωs[i + 1]
    frac = tr - i + 1
    return Ω1 * (1 - frac) + Ω2 * frac
end

function combine_candidates(cands, weights)
    ts = [0.0]
    nts = 1
    for c in cands
        n = length(c.Ωs) - 1
        new_nts = nts + n
        resize!(ts, new_nts)
        ts[nts + 1:new_nts] .= (1:n) .* c.τ
        nts = new_nts
    end
    sort!(ts)
    vals = zeros(nts)
    for (ci, c) in enumerate(cands)
        w = weights[ci]
        for (i, t) in enumerate(ts)
            vals[i] += w * eval_cand(c, t)
        end
    end
    return ts, vals
end

const NIons = 13
const tgt_scale = 0.1
const tgt = zeros(13, 13)
for i in 2:NIons - 1
    for j in 2:NIons - 1
        if abs(i - j) <= 2 && i != j
            tgt[i, j] = 0.1 * tgt_scale
        end
    end
end
# for i in 1:NIons - 3
#     tgt[i + 3, i] = tgt[i, i + 3] = tgt_scale * (((i % 4) in (1, 2)) ? 1 : 2)
# end
# for i in 1:8:NIons
#     tgt[i + 2, i + 3] = tgt[i + 3, i + 2] = tgt[i, i + 1] = tgt[i + 1, i] = tgt_scale * 0.5
# end

const kernel = AreaKernel(candidates, sysparams)
const tracker = get_tracker(kernel, 0.2)
const opt = get_opt(kernel, tgt, tracker)

val, args = opt_rep(opt, tracker, 10)

using NaCsPlot
using PyPlot

for n in 1:11
    sub_args = @view args[kernel.ncands * n + 1:kernel.ncands * n + kernel.ncands]
    ts, vs = combine_candidates(candidates, sub_args)
    @show maximum(abs.(vs))
    plot(ts, vs .+ (n - 1) * 0.5)
end
# yticks((0:10) .* 2)
# grid()
# show()

#

include("am_shared.jl")
include("ind-modulation.jl")

using LinearAlgebra
using NLopt
using MSSim: Optimizers as Opts
using GoldGates

smooth_spec(nseg; width, seg_samples) =
    Seq.ModSpec{nseg * seg_samples}(amp=Seq.AmpSpec(cb=get_smooth_am_cbs(nseg, width), sym=false))

function constraint_matrix(spec::Seq.ModSpec{_nseg}, τ, δ, ωms; seg_samples) where _nseg
    param_buff = zeros(Seq.nparams(spec))
    @inbounds param_buff[spec.τ] = τ / seg_samples
    @inbounds for i in spec.ωs
        param_buff[i] = δ
    end
    @inbounds param_buff[spec.Ωs[1]] = 1

    raw_param_buff = zeros(_nseg * 5)
    Seq.transform_argument(spec, raw_param_buff, param_buff)

    nmodes = length(ωms)
    namps = length(spec.Ωs)
    @assert namps > 4 * nmodes
    res = Matrix{Float64}(undef, 4 * nmodes, namps)

    mask = SS.ValueMask(true, true, false, false, true, false)
    buf = SL.ComputeBuffer{_nseg,Float64}(Val(mask), Val(zero(SS.ValueMask)))
    kern = SL.Kernel(buf, Val(zero(SL.ParamGradMask)))

    for (ωi, ωm) in enumerate(ωms)
        SL.eval_with_mode!(kern, raw_param_buff, ωm)
        ω = δ - ωm
        val0 = kern.result.val
        for Ωi in 1:namps
            t0 = (Ωi - 1) * τ
            dis = val0.dis * cis(ω * t0)
            disδ = val0.disδ * cis(ω * t0) + im * t0 * dis
            res[ωi * 4 - 3, Ωi] = real(dis)
            res[ωi * 4 - 2, Ωi] = imag(dis)
            res[ωi * 4 - 1, Ωi] = real(disδ)
            res[ωi * 4 - 0, Ωi] = imag(disδ)
        end
    end
    return res
end

function area2_matrix(spec::Seq.ModSpec{_nseg}, τ, δ, ωm; width, seg_samples, areas=nothing) where _nseg
    namps = length(spec.Ωs)
    if areas === nothing
        areas = Matrix{Float64}(undef, namps, namps)
    else
        @assert size(areas) == (namps, namps)
    end
    @inbounds for i in 1:namps
        a1 = spec.amps[i]
        for j in 1:namps
            a2 = spec.amps[j]
            areas[j, i] = compute_area2_am(τ / seg_samples, τ / seg_samples,
                                           δ, ωm, a1, a2)
        end
    end
    @inbounds for i in 1:namps
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

struct AreaKernel{NSeg,RampWidth,SegSamples}
    nions::Int
    ncands::Int
    waveforms::Matrix{Float64}
    all_areas::Vector{Matrix{Float64}}
    all_couple::Vector{Matrix{Float64}}
    cs′::Vector{Matrix{Float64}}
    function AreaKernel{NSeg,RampWidth,SegSamples}(sysparams, τ, δ) where {NSeg,RampWidth,SegSamples}
        ωs = 2π .* sysparams.modes.radial1
        nions = length(ωs)
        _nseg = NSeg * SegSamples
        spec = smooth_spec(NSeg; width=RampWidth, seg_samples=SegSamples)
        C = constraint_matrix(spec, τ, δ, ωs; seg_samples=SegSamples)
        candidates = nullspace(C)
        ncands = size(candidates, 2)
        namps = length(spec.Ωs)

        orig_basis = Matrix{Float64}(undef, _nseg + 1, namps)
        @inbounds for i in 1:namps
            orig_basis[:, i] .= spec.amps[i]
        end
        waveforms = orig_basis * candidates

        all_areas = Matrix{Float64}[]
        areas_buff = Matrix{Float64}(undef, namps, namps)
        for ωm in ωs
            area2_matrix(spec, τ, δ, ωm; width=RampWidth, seg_samples=SegSamples,
                         areas=areas_buff)
            push!(all_areas, candidates' * areas_buff * candidates)
        end

        return new{NSeg,RampWidth,SegSamples}(
            nions,
            ncands,
            waveforms,
            all_areas,
            [couple_matrix(sysparams, modei, nions) for modei in 1:nions],
            [Matrix{Float64}(undef, ncands, nions) for modei in 1:nions])
    end
end

get_waveform(k::AreaKernel, weights) =k.waveforms * weights

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

params_file = "072125_goldparams_13ions.json"
sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
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

const kernel = AreaKernel{200, 2, 10}(sysparams, 1,
                                      2π * (sysparams.modes.radial1[1] - 0.01))
const tracker = get_tracker(kernel, 0.2)
const opt = get_opt(kernel, tgt, tracker)

val, args = opt_rep(opt, tracker, 10)

using PyPlot

for n in 1:11
    sub_args = @view args[kernel.ncands * n + 1:kernel.ncands * n + kernel.ncands]
    vs = get_waveform(kernel, sub_args)
    @show maximum(abs.(vs))
    plot(vs .+ (n - 1) * 0.5)
end

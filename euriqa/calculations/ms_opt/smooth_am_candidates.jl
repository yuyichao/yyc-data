#

include("am_shared.jl")
include("ind-modulation.jl")

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

function area2_matrix(spec::Seq.ModSpec{_nseg}, τ, δ, ωm; width, seg_samples) where _nseg
    namps = length(spec.Ωs)
    areas = Matrix{Float64}(undef, namps, namps)
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

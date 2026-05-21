#!/usr/bin/julia

using MSSim: Utils as U, SegSeq as SS, SymLinear as SL

const rf_vv = U.cos_f1
const rf_vd = U.sin_f1
const rf_dv = U.TrigRatio{true,3,(),(1,),(-1,)}()
const rf_dd = U.cos_f3

# const if_vv = U.TrigRatio{true,2,(1,),(-1,),()}()
# const if_vd = U.TrigRatio{false,3,(-1,1//2),(),(1,)}()
# const if_dv = U.TrigRatio{false,3,(1,1//2),(-1,),(-1,)}()
# const if_dd = U.sin_f3

@inline function enclosed_area2_kernel(o1, o1′, o2, o2′, d, s, c, V)
    c_vv = V.C1
    c_vd = V.S1
    c_dv = @inline rf_dv(d, s, c)
    c_dd = V.C3

    vv = o1 * o2
    vd = o1 * o2′
    dv = o1′ * o2
    dd = o1′ * o2′

    return muladd(c_vv,  vv, muladd(c_vd, vd, muladd(c_dv, dv, c_dd * dd)))
end

struct SegData2{T}
    τ::T
    dis1::Complex{T}
    dis2::Complex{T}
    area::T
end

@inline function compute_values2(τ::_T, Ω1, Ω1′, Ω2, Ω2′, φ, δ)
    T = float(_T)
    CT = Complex{T}
    SDV = SegData2(T)

    @inline begin
        d = δ * τ
        o1 = Ω1 * τ
        o1′ = Ω1′ * τ^2
        o2 = Ω2 * τ
        o2′ = Ω2′ * τ^2
        s, c = U.fast_sincos(d)
        sφ, cφ = U.fast_sincos(φ)
        phase0 = complex(cφ, sφ)
        phase0_τ = phase0 * τ
        V = SL.@gen_trig_ratios(d, s, c, sin_c1, sin_c2, sin_f1, cos_f1, cos_f3)

        dis1 = U.mul(phase0, SL.displacement_kernel(o1, o1′, d, s, c, V))
        dis2 = U.mul(phase0, SL.displacement_kernel(o2, o2′, d, s, c, V))
        area = enclosed_area2_kernel(o1, o1′, o2, o2′, d, s, c, V)
        res = SDV(τ, dis1, dis2, area)
    end
    return res
end

function compute_area2_single_mode(segments::AbstractVector{SegData2{T}}) where {T}
    nseg = length(segments)

    p_τ = zero(T)
    p_dis1 = complex(zero(T))
    p_dis2 = complex(zero(T))
    p_area = zero(T)
    @inbounds for i in 1:nseg
        seg = segments[i]

        np_τ = p_τ + seg.τ
        np_dis1 = p_dis1 + seg.dis1
        np_dis2 = p_dis2 + seg.dis2
        np_area = p_area + muladd(real(p_dis2), imag(seg.dis1),
                                  muladd(-imag(p_dis2), real(seg.dis1), seg.area))

        p_τ = np_τ
        p_dis1 = np_dis1
        p_dis2 = np_dis2
        p_area = np_area
    end

    return SegData2{T}(p_τ, p_dis1, p_dis2, p_area)
end

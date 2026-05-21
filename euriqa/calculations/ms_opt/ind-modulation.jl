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

Base.zero(::Type{SegData2{T}}) where T = SegData2{T}(zero(T), zero(Complex{T}),
                                                     zero(Complex{T}), zero(T))

@inline function compute_values2(τ::_T, Ω1, Ω1′, Ω2, Ω2′, φ, δ) where _T
    T = float(_T)
    CT = Complex{T}
    SDV = SegData2{T}

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
        V = SL.SegInt.@gen_trig_ratios(d, s, c, sin_c1, sin_c2, sin_f1, cos_f1, cos_f3)

        dis1 = U.mul(phase0, SL.SegInt.displacement_kernel(o1, o1′, d, s, c, V))
        dis2 = U.mul(phase0, SL.SegInt.displacement_kernel(o2, o2′, d, s, c, V))
        area = enclosed_area2_kernel(o1, o1′, o2, o2′, d, s, c, V)
        res = SDV(τ, dis1, dis2, area)
    end
    return res
end

@inline function add_segment(cum::SegData2{T}, seg::SegData2{T}) where T
    return SegData2{T}(cum.τ + seg.τ,
                       cum.dis1 + seg.dis1,
                       cum.dis2 + seg.dis2,
                       cum.area + muladd(real(cum.dis2), imag(seg.dis1),
                                         muladd(-imag(cum.dis2), real(seg.dis1),
                                                seg.area)))
end

@inline function compute_area2(args, ωm)
    T = eltype(args)
    ωm = T(ωm)

    narg = length(args)
    nseg = narg ÷ 7
    @assert narg == nseg * 7

    φm = zero(T)
    res = zero(SegData2{T})

    @inbounds for i in 1:nseg
        τ = args[i * 7 - 6]
        Ω1 = args[i * 7 - 5]
        Ω1′ = args[i * 7 - 4]
        Ω2 = args[i * 7 - 3]
        Ω2′ = args[i * 7 - 2]
        φ = args[i * 7 - 1] - φm
        δ = args[i * 7] - ωm

        φm = muladd(ωm, τ, φm)

        res = add_segment(res, compute_values2(τ, Ω1, Ω1′, Ω2, Ω2′, φ, δ))
    end
    return res.area
end

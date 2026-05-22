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

struct AreaModeState{T}
    φm::T
    res::SegData2{T}
end
Base.zero(::Type{AreaModeState{T}}) where T = AreaModeState{T}(zero(T), zero(SegData2{T}))

@inline function add_segment(state::AreaModeState{T}, τ, Ω1, Ω1′, Ω2, Ω2′, φ, δ, ωm) where T
    φm = state.φm
    res = state.res

    φ -= φm
    δ -= ωm

    φm = muladd(ωm, τ, φm)

    res = add_segment(res, compute_values2(τ, Ω1, Ω1′, Ω2, Ω2′, φ, δ))

    return AreaModeState{T}(φm, res)
end

mutable struct SeqStatus{T}
    i::Int # pointing to end/next start
    const nΩ::Int
    tend::T
    Ω::T
    Ωe::T
    Ω′::T

    function SeqStatus{T}(dτ, Ωs) where T
        Ω = Ωs[1]
        Ωe = Ωs[2]
        return new{T}(2, length(Ωs), dτ, Ω, Ωe, (Ωe - Ω) / dτ)
    end
end

@inline function next_step(s::SeqStatus{T}, dτ, Ωs, tstart, tend) where T
    if s.i > s.nΩ
        # Sequence already finished
        return
    end
    if tend < s.tend
        # Step on-going
        s.Ω = muladd(s.Ω′, tend - tstart, s.Ω)
        return
    end
    if s.i == s.nΩ
        # Sequence finishing
        s.tend = Inf
        s.Ω = 0.0
        s.Ωe = 0.0
        s.Ω′ = 0.0
    else
        s.tend = s.i * dτ
        s.Ω = s.Ωe
        s.Ωe = Ωs[s.i + 1]
        s.Ω′ = (s.Ωe - s.Ω) / dτ
    end
    s.i += 1
    return
end

function compute_area2_am(dτ1, dτ2, ω, ωm, Ω1s, Ω2s)
    state = zero(AreaModeState{Float64})
    s1 = SeqStatus{Float64}(dτ1, Ω1s)
    s2 = SeqStatus{Float64}(dτ2, Ω2s)

    tstart = 0.0
    tend = min(s1.tend, s2.tend)

    while s1.i <= s1.nΩ || s2.i <= s2.nΩ
        τ = tend - tstart
        state = add_segment(state, τ, s1.Ω, s1.Ω′, s2.Ω, s2.Ω′, ω * tstart, ω, ωm)
        next_step(s1, dτ1, Ω1s, tstart, tend)
        next_step(s2, dτ2, Ω2s, tstart, tend)
        tstart = tend
        tend = min(s1.tend, s2.tend)
    end

    return state.res.area
end

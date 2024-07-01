#!/usr/bin/julia

using QuantumOptics
using LinearAlgebra

struct LabBasisData{B,K,O} <: Function
    Ω10_2::Float64
    dΩ1_2::Float64
    Ω20_2::Float64
    dΩ2_2::Float64
    δ10::Float64
    dδ1_2::Float64
    δ20::Float64
    dδ2_2::Float64
    tlen::Float64

    basis::B
    ψ0::K
    op::O
    function LabBasisData(Ω10, Ω11, Ω20, Ω21, δ10, δ11, δ20, δ21, tlen)
        basis = GenericBasis(3)
        op_data = zeros(ComplexF64, 3, 3)
        op = Operator(basis, op_data)
        ψ0 = basisstate(basis, 1)
        Ω10_2 = Ω10 / 2
        Ω20_2 = Ω20 / 2
        dΩ1_2 = (Ω11 - Ω10) / tlen / 2
        dΩ2_2 = (Ω21 - Ω20) / tlen / 2
        dδ1_2 = (δ11 - δ10) / tlen / 2
        dδ2_2 = (δ21 - δ20) / tlen / 2
        return new{typeof(basis),typeof(ψ0),typeof(op)}(Ω10_2, dΩ1_2, Ω20_2, dΩ2_2,
                                                         δ10, dδ1_2, δ20, dδ2_2,
                                                         tlen, basis, ψ0, op)
    end
end

function (data::LabBasisData)(t, ψ)
    op_data = data.op.data
    # δi = δi0 + (δi1 - δi0) * t / tlen
    # φi = ∫δidt = δi0 * t + (δi1 - δi0) * t^2 / 2 / tlen
    φ1 = t * muladd(data.dδ1_2, t, data.δ10)
    φ2 = t * muladd(data.dδ2_2, t, data.δ20)
    Ω1_2 = muladd(data.dΩ1_2, t, data.Ω10_2)
    Ω2_2 = muladd(data.dΩ2_2, t, data.Ω20_2)
    M1 = Ω1_2 * cis(φ1)
    M2 = Ω2_2 * cis(φ2)
    op_data[1, 2] = M1
    op_data[2, 1] = conj(M1)
    op_data[2, 3] = M2
    op_data[3, 2] = conj(M2)
    return data.op
end

struct LabBasisCalc{LD}
    data::LD
    function LabBasisCalc(args...)
        ld = LabBasisData(args...)
        return new{typeof(ld)}(ld)
    end
end

function evolve(calc::LabBasisCalc, npoints=1001; kws...)
    calc_data = calc.data
    ts, ψs = timeevolution.schroedinger_dynamic(range(0, calc_data.tlen, npoints),
                                                 calc_data.ψ0, calc_data; kws...)
    return ts, [ψ.data[i] for i in 1:3, ψ in ψs]
end

struct RotBasisData{B,K,O} <: Function
    Ω10_2::Float64
    dΩ1_2::Float64
    Ω20_2::Float64
    dΩ2_2::Float64
    δ10::Float64
    dδ1_2::Float64
    δ20::Float64
    dδ2_2::Float64

    tlen::Float64
    basis::B
    ψ0::K
    op::O
    function RotBasisData(Ω10, Ω11, Ω20, Ω21, δ10, δ11, δ20, δ21, tlen)
        basis = GenericBasis(3)
        op_data = zeros(ComplexF64, 3, 3)
        op = Operator(basis, op_data)
        ψ0 = basisstate(basis, 1)
        Ω10_2 = Ω10 / 2
        Ω20_2 = Ω20 / 2
        dΩ1_2 = (Ω11 - Ω10) / tlen / 2
        dΩ2_2 = (Ω21 - Ω20) / tlen / 2
        dδ1_2 = (δ11 - δ10) / tlen / 2
        dδ2_2 = (δ21 - δ20) / tlen / 2
        return new{typeof(basis),typeof(ψ0),typeof(op)}(Ω10_2, dΩ1_2, Ω20_2, dΩ2_2,
                                                         δ10, dδ1_2, δ20, dδ2_2,
                                                         tlen, basis, ψ0, op)
    end
end

function (data::RotBasisData)(t, ψ)
    op_data = data.op.data
    δ1 = muladd(data.dδ1_2 * 2, t, data.δ10)
    δ2 = muladd(data.dδ2_2 * 2, t, data.δ20)
    Ω1_2 = muladd(data.dΩ1_2, t, data.Ω10_2)
    Ω2_2 = muladd(data.dΩ2_2, t, data.Ω20_2)
    op_data[1, 1] = δ1
    op_data[3, 3] = -δ2
    op_data[1, 2] = Ω1_2
    op_data[2, 1] = Ω1_2
    op_data[2, 3] = Ω2_2
    op_data[3, 2] = Ω2_2
    return data.op
end

struct RotBasisCalc{LD}
    data::LD
    function RotBasisCalc(args...)
        ld = RotBasisData(args...)
        return new{typeof(ld)}(ld)
    end
end

function evolve(calc::RotBasisCalc, npoints=1001; kws...)
    calc_data = calc.data
    ts, ψs = timeevolution.schroedinger_dynamic(range(0, calc_data.tlen, npoints),
                                                 calc_data.ψ0, calc_data; kws...)
    data = [ψ.data[i] for i in 1:3, ψ in ψs]
    δ10 = calc_data.δ10
    δ20 = calc_data.δ20
    dδ1_2 = calc_data.dδ1_2
    dδ2_2 = calc_data.dδ2_2
    tlen = calc_data.tlen
    for (i, t) in enumerate(ts)
        φ1 = t * muladd(dδ1_2, t, δ10)
        φ2 = t * muladd(dδ2_2, t, δ20)
        p1 = cis(φ1)
        p2 = cis(φ2)
        data[1, i] *= p1
        data[3, i] *= conj(p2)
    end
    return ts, data
end

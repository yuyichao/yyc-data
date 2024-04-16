#!/usr/bin/julia

using QuantumOptics

struct LabBasisData{B,K,O} <: Function
    Ω0_2::Float64
    dΩ_2::Float64
    δ0::Float64
    dδ_2::Float64
    tlen::Float64

    basis::B
    ψ0::K
    op::O
    function LabBasisData(Ω0, Ω1, δ0, δ1, tlen)
        basis = SpinBasis(1//2)
        op_data = zeros(ComplexF64, 2, 2)
        op = Operator(basis, op_data)
        ψ0 = spindown(basis)
        Ω0_2 = Ω0 / 2
        dΩ_2 = (Ω1 - Ω0) / tlen / 2
        dδ_2 = (δ1 - δ0) / tlen / 2
        return new{typeof(basis),typeof(ψ0),typeof(op)}(Ω0_2, dΩ_2, δ0, dδ_2,
                                                         tlen, basis, ψ0, op)
    end
end

function (data::LabBasisData)(t, ψ)
    op_data = data.op.data
    # δ = δ0 + (δ1 - δ0) * t / tlen
    # φ = ∫δdt = δ0 * t + (δ1 - δ0) * t^2 / 2 / tlen
    φ = t * muladd(data.dδ_2, t, data.δ0)
    Ω_2 = muladd(data.dΩ_2, t, data.Ω0_2)
    M = Ω_2 * cis(φ)
    op_data[1, 2] = M
    op_data[2, 1] = conj(M)
    return data.op
end

struct LabBasisCalc{LD}
    data::LD
    function LabBasisCalc(args...)
        ld = LabBasisData(args...)
        return new{typeof(ld)}(ld)
    end
end

function evolve(calc::LabBasisCalc, npoints=1001)
    calc_data = calc.data
    ts, ψs = timeevolution.schroedinger_dynamic(range(0, calc_data.tlen, npoints),
                                                 calc_data.ψ0, calc_data)
    return ts, [ψ.data[i] for i in 1:2, ψ in ψs]
end

struct RotBasisData{B,K,O} <: Function
    Ω0_2::Float64
    dΩ_2::Float64
    δ0_2::Float64
    dδ_2::Float64
    tlen::Float64
    basis::B
    ψ0::K
    op::O
    function RotBasisData(Ω0, Ω1, δ0, δ1, tlen)
        basis = SpinBasis(1//2)
        op_data = zeros(ComplexF64, 2, 2)
        op = Operator(basis, op_data)
        ψ0 = spindown(basis)
        Ω0_2 = Ω0 / 2
        dΩ_2 = (Ω1 - Ω0) / tlen / 2
        dδ_2 = (δ1 - δ0) / tlen / 2
        return new{typeof(basis),typeof(ψ0),typeof(op)}(Ω0_2, dΩ_2, δ0 / 2, dδ_2, tlen,
                                                         basis, ψ0, op)
    end
end

function (data::RotBasisData)(t, ψ)
    op_data = data.op.data
    δ_2 = muladd(data.dδ_2, t, data.δ0_2)
    Ω_2 = muladd(data.dΩ_2, t, data.Ω0_2)
    op_data[1, 1] = δ_2
    op_data[2, 2] = -δ_2
    op_data[1, 2] = Ω_2
    op_data[2, 1] = Ω_2
    return data.op
end

struct RotBasisCalc{LD}
    data::LD
    function RotBasisCalc(args...)
        ld = RotBasisData(args...)
        return new{typeof(ld)}(ld)
    end
end

function evolve(calc::RotBasisCalc, npoints=1001)
    calc_data = calc.data
    ts, ψs = timeevolution.schroedinger_dynamic(range(0, calc_data.tlen, npoints),
                                                 calc_data.ψ0, calc_data)
    data = [ψ.data[i] for i in 1:2, ψ in ψs]
    δ0_2 = calc_data.δ0_2
    dδ_4 = calc_data.dδ_2 / 2
    tlen = calc_data.tlen
    for (i, t) in enumerate(ts)
        φ_2 = t * muladd(dδ_4, t, δ0_2)
        p = cis(φ_2)
        data[1, i] *= p
        data[2, i] *= conj(p)
    end
    return ts, data
end

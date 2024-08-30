#!/usr/bin/julia

using QuantumOptics
using LinearAlgebra
using CGcoefficient

function g_sum(J, J1, g1, J2, g2)
    return (g1 * (J * (J + 1) + J1 * (J1 + 1) - J2 * (J2 + 1)) + g2 * (J * (J + 1) + J2 * (J2 + 1) - J1 * (J1 + 1))) / (2 * J * (J + 1))
end

const g_s = 2.0023
const g_l = 1.0

const g_S1_2 = g_s
const g_P1_2 = g_sum(1/2, 1, g_l, 1/2, g_s)
const g_D3_2 = g_sum(3/2, 2, g_l, 1/2, g_s)
const g_F7_2 = g_sum(7/2, 3, g_l, 1/2, g_s)
const g_B1_2 = g_sum(1/2, 3/2, g_sum(3/2, 7/2, g_F7_2, 2, g_l), 1, g_s)

const μB = 2π * 1.4e6

const ΓP_S = 1.23e8
const ΓP_D = 6.19e5
const ΓB_S = 2.61e7
const ΓB_D = 4.75e5
# const ΓD_S = 18.98

const HF_S1_2 = 2π * 12.643e9
const HF_P1_2 = 2π * 2.105e9
const HF_D3_2 = 2π * 0.86e9
const HF_B = -2π * 2.2095e9

function get_hf_hamiltonian(gJ, Fl, hfl, hfh)
    dJ = 2 * Fl + 1
    J = dJ / 2
    dI = 1

    Fh = Fl + 1
    nFl = 2 * Fl + 1
    nFh = 2 * Fh + 1

    bl = SpinBasis(Fl)
    bh = SpinBasis(Fh)
    b = bl ⊕ bh

    op_data = zeros(Float64, nFl + nFh, nFl + nFh)

    # stretched states
    op_data[nFl + 1, nFl + 1] = -J * gJ + hfh
    op_data[nFl + nFh, nFl + nFh] = J * gJ + hfh

    for mF in -Fl:Fl
        il = mF + Fl + 1
        ih = mF + Fl + nFl + 2
        dmj1 = 2 * mF - 1
        dmj2 = 2 * mF + 1
        mj1 = dmj1 / 2
        mj2 = dmj2 / 2

        cg1l = fCG(dJ, dI, Fl * 2, dmj1, 1, mF * 2)
        cg2l = fCG(dJ, dI, Fl * 2, dmj2, -1, mF * 2)
        cg1h = fCG(dJ, dI, Fh * 2, dmj1, 1, mF * 2)
        cg2h = fCG(dJ, dI, Fh * 2, dmj2, -1, mF * 2)

        op_data[il, il] = gJ * (mj1 * cg1l^2 + mj2 * cg2l^2) + hfl
        op_data[il, ih] = op_data[ih, il] = gJ * (mj1 * cg1l * cg1h + mj2 * cg2l * cg2h)
        op_data[ih, ih] = gJ * (mj1 * cg1h^2 + mj2 * cg2h^2) + hfh
    end

    return Operator(b, op_data)
end

# State order,
# S1/2 (F=0, F=1), P1/2 (F=0, F=1), D3/2 (F=1, F=2), [3/2]1/2 (F=0, F=1)

struct Yb171SysData{Basis,O} <: Function
    basis::Basis
    op::O
    h0::H
    function Yb171SysData(B)
        h_S1 = get_hf_hamiltonian(g_S1_2 * μB * B, 0, -HF_S1_2, 0)
        h_P1 = get_hf_hamiltonian(g_P1_2 * μB * B, 0, 0, HF_P1_2)
        h_D3 = get_hf_hamiltonian(g_D3_2 * μB * B, 1, 0, HF_D3_2)
        h_B = get_hf_hamiltonian(g_B1_2 * μB * B, 0, 0, HF_B)

        h0 = h_S1 ⊕ h_P1 ⊕ h_D3 ⊕ h_B

        # Γσ⁺
        # Γσ⁻
        # Γπ
    end
end

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

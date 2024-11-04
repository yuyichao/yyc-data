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

function dipole_branch(q, dJ, dJ′, dI, dF, dF′, dmF, dmF′)
    return ((-1)^(dF′ + (dmF + dJ + dI) ÷ 2) *
        sqrt((dJ′ + 1) * (dF + 1) * (dF′ + 1)) *
        f6j(dJ, dJ′, 2, dF′, dF, dI) * f3j(dF′, 2, dF, dmF′, 2 * q, dmF))
end

function fill_J!(op, Γ, q, dJ, dJ′, dI, idxFl, idxFl′)
    sqrtΓ = sqrt(Γ)
    dFs = abs(dJ - dI):2:(dJ + dI)
    dF′s = abs(dJ′ - dI):2:(dJ′ + dI)
    for (i, dF) in enumerate(dFs)
        idxF = idxFl + i - 1
        basisF = SpinBasis(dF//2)
        for (i′, dF′) in enumerate(dF′s)
            idxF′ = idxFl′ + i′ - 1
            basisF′ = SpinBasis(dF′//2)
            subop = Operator(basisF, basisF′,
                             [sqrtΓ * dipole_branch(q, dJ, dJ′, dI, dF, dF′, dmF, dmF′)
                              for dmF in -dF:2:dF, dmF′ in -dF′:2:dF′])
            setblock!(op, subop, idxF, idxF′)
        end
    end
    return op
end

function fill_dipole!(op⁺, op⁻, q, dJ, dJ′, dI, idxFl, idxFl′)
    dFs = abs(dJ - dI):2:(dJ + dI)
    dF′s = abs(dJ′ - dI):2:(dJ′ + dI)
    for (i, dF) in enumerate(dFs)
        idxF = idxFl + i - 1
        basisF = SpinBasis(dF//2)
        for (i′, dF′) in enumerate(dF′s)
            idxF′ = idxFl′ + i′ - 1
            basisF′ = SpinBasis(dF′//2)
            subop = Operator(basisF, basisF′,
                             [dipole_branch(q, dJ, dJ′, dI, dF, dF′, dmF, dmF′)
                              for dmF in -dF:2:dF, dmF′ in -dF′:2:dF′])
            setblock!(op⁻, subop, idxF, idxF′)
            setblock!(op⁺, subop', idxF′, idxF)
        end
    end
    return op⁺, op⁻
end

# State order,
# S1/2 (F=0, F=1), P1/2 (F=0, F=1), D3/2 (F=1, F=2), [3/2]1/2 (F=0, F=1)
#         1,   2,          3,   4,          5,   6,              7,   8

struct Yb171Sys{Basis,O}
    basis::Basis
    op::O
    h0::O

    J::Vector{O}
    Jdagger::Vector{O}

    dP_Sσ⁺⁺::O
    dP_Sσ⁺⁻::O
    dP_Sσ⁻⁺::O
    dP_Sσ⁻⁻::O
    dP_Sπ⁺::O
    dP_Sπ⁻::O

    dP_Dσ⁺⁺::O
    dP_Dσ⁺⁻::O
    dP_Dσ⁻⁺::O
    dP_Dσ⁻⁻::O
    dP_Dπ⁺::O
    dP_Dπ⁻::O

    dB_Sσ⁺⁺::O
    dB_Sσ⁺⁻::O
    dB_Sσ⁻⁺::O
    dB_Sσ⁻⁻::O
    dB_Sπ⁺::O
    dB_Sπ⁻::O

    dB_Dσ⁺⁺::O
    dB_Dσ⁺⁻::O
    dB_Dσ⁻⁺::O
    dB_Dσ⁻⁻::O
    dB_Dπ⁺::O
    dB_Dπ⁻::O
    function Yb171Sys(B)
        h_S1 = get_hf_hamiltonian(g_S1_2 * μB * B, 0, -HF_S1_2, 0)
        h_P1 = get_hf_hamiltonian(g_P1_2 * μB * B, 0, 0, HF_P1_2)
        h_D3 = get_hf_hamiltonian(g_D3_2 * μB * B, 1, 0, HF_D3_2)
        h_B = get_hf_hamiltonian(g_B1_2 * μB * B, 0, 0, HF_B)

        h0 = h_S1 ⊕ h_P1 ⊕ h_D3 ⊕ h_B

        basis = h0.basis_l

        J = [fill_J!(zero(h0), ΓP_S, 1, 1, 1, 1, 1, 3),
             fill_J!(zero(h0), ΓP_S, -1, 1, 1, 1, 1, 3),
             fill_J!(zero(h0), ΓP_S, 0, 1, 1, 1, 1, 3),
             fill_J!(zero(h0), ΓP_D, 1, 3, 1, 1, 5, 3),
             fill_J!(zero(h0), ΓP_D, -1, 3, 1, 1, 5, 3),
             fill_J!(zero(h0), ΓP_D, 0, 3, 1, 1, 5, 3),
             fill_J!(zero(h0), ΓB_S, 1, 1, 1, 1, 1, 7),
             fill_J!(zero(h0), ΓB_S, -1, 1, 1, 1, 1, 7),
             fill_J!(zero(h0), ΓB_S, 0, 1, 1, 1, 1, 7),
             fill_J!(zero(h0), ΓB_D, 1, 3, 1, 1, 5, 7),
             fill_J!(zero(h0), ΓB_D, -1, 3, 1, 1, 5, 7),
             fill_J!(zero(h0), ΓB_D, 0, 3, 1, 1, 5, 7)]
        Jdagger = dagger.(J)

        dP_Sσ⁺⁺, dP_Sσ⁺⁻ = fill_dipole!(zero(h0), zero(h0), 1, 1, 1, 1, 1, 3)
        dP_Sσ⁻⁺, dP_Sσ⁻⁻ = fill_dipole!(zero(h0), zero(h0), -1, 1, 1, 1, 1, 3)
        dP_Sπ⁺, dP_Sπ⁻ = fill_dipole!(zero(h0), zero(h0), 0, 1, 1, 1, 1, 3)

        dP_Dσ⁺⁺, dP_Dσ⁺⁻ = fill_dipole!(zero(h0), zero(h0), 1, 3, 1, 1, 5, 3)
        dP_Dσ⁻⁺, dP_Dσ⁻⁻ = fill_dipole!(zero(h0), zero(h0), -1, 3, 1, 1, 5, 3)
        dP_Dπ⁺, dP_Dπ⁻ = fill_dipole!(zero(h0), zero(h0), 0, 3, 1, 1, 5, 3)

        dB_Sσ⁺⁺, dB_Sσ⁺⁻ = fill_dipole!(zero(h0), zero(h0), 1, 1, 1, 1, 1, 7)
        dB_Sσ⁻⁺, dB_Sσ⁻⁻ = fill_dipole!(zero(h0), zero(h0), -1, 1, 1, 1, 1, 7)
        dB_Sπ⁺, dB_Sπ⁻ = fill_dipole!(zero(h0), zero(h0), 0, 1, 1, 1, 1, 7)

        dB_Dσ⁺⁺, dB_Dσ⁺⁻ = fill_dipole!(zero(h0), zero(h0), 1, 3, 1, 1, 5, 7)
        dB_Dσ⁻⁺, dB_Dσ⁻⁻ = fill_dipole!(zero(h0), zero(h0), -1, 3, 1, 1, 5, 7)
        dB_Dπ⁺, dB_Dπ⁻ = fill_dipole!(zero(h0), zero(h0), 0, 3, 1, 1, 5, 7)

        return new{typeof(basis),typeof(h0)}(
            basis, zero(h0), h0, J, Jdagger,
            dP_Sσ⁺⁺, dP_Sσ⁺⁻, dP_Sσ⁻⁺, dP_Sσ⁻⁻, dP_Sπ⁺, dP_Sπ⁻,
            dP_Dσ⁺⁺, dP_Dσ⁺⁻, dP_Dσ⁻⁺, dP_Dσ⁻⁻, dP_Dπ⁺, dP_Dπ⁻,
            dB_Sσ⁺⁺, dB_Sσ⁺⁻, dB_Sσ⁻⁺, dB_Sσ⁻⁻, dB_Sπ⁺, dB_Sπ⁻,
            dB_Dσ⁺⁺, dB_Dσ⁺⁻, dB_Dσ⁻⁺, dB_Dσ⁻⁻, dB_Dπ⁺, dB_Dπ⁻
        )
    end
end

function get_ψ(sys::Yb171Sys, name::Symbol, F, mF)
    ψ = Ket(sys.basis)
    @assert -F <= mF <= F
    if name === :S
        @assert 0 <= F <= 1
        ψ.data[F * 2 + mF + 1] = 1
    elseif name === :P
        @assert 0 <= F <= 1
        ψ.data[F * 2 + mF + 5] = 1
    elseif name === :D
        @assert 1 <= F <= 2
        ψ.data[F * 4 + mF + 6] = 1
    elseif name === :B
        @assert 0 <= F <= 1
        ψ.data[F * 2 + mF + 17] = 1
    else
        @assert false
    end
    return ψ
end

function get_ρ(sys::Yb171Sys, name::Symbol, F, mF)
    ψ = get_ψ(sys, name, F, mF)
    return ψ ⊗ ψ'
end

function evolve(update!, sys::Yb171Sys, ρ0, tlen, npoints=1001; kws...)
    cb = let update! = update!, sys = sys
        function (t, rho)
            update!(sys, t)
            return (sys.op, sys.J, sys.Jdagger)
        end
    end
    return timeevolution.master_dynamic(range(0, tlen, npoints), ρ0, cb; kws...)
end

function no_drive!(sys::Yb171Sys, t)
    sys.op .= sys.h0
end

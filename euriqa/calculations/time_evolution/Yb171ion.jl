#!/usr/bin/julia

using QuantumOptics
using LinearAlgebra
using CGcoefficient
using SparseArrays: sparse, spzeros, SparseMatrixCSC

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

    Idx = Int[]
    Idx′ = Int[]
    V = ComplexF64[]

    function add_element!(i, j, v)
        push!(Idx, i)
        push!(Idx′, j)
        push!(V, v)
        return v
    end

    # stretched states
    add_element!(nFl + 1, nFl + 1, -J * gJ + hfh)
    add_element!(nFl + nFh, nFl + nFh, J * gJ + hfh)

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

        add_element!(il, il, gJ * (mj1 * cg1l^2 + mj2 * cg2l^2) + hfl)
        add_element!(ih, il,
                     add_element!(il, ih, gJ * (mj1 * cg1l * cg1h + mj2 * cg2l * cg2h)))
        add_element!(ih, ih, gJ * (mj1 * cg1h^2 + mj2 * cg2h^2) + hfh)
    end

    return Operator(b, sparse(Idx, Idx′, V, nFl + nFh, nFl + nFh))
end

function dipole_branch(q, dJ, dJ′, dI, dF, dF′, dmF, dmF′)
    return ((-1)^(dF′ + (dmF + dJ + dI) ÷ 2) *
        sqrt((dJ′ + 1) * (dF + 1) * (dF′ + 1)) *
        f6j(dJ, dJ′, 2, dF′, dF, dI) * f3j(dF′, 2, dF, dmF′, 2 * q, -dmF))
end

function add_J!(Js, basis, offsets, Γ, dJ, dJ′, dI, idxFl, idxFl′)
    sqrtΓ = sqrt(Γ)
    dFs = abs(dJ - dI):2:(dJ + dI)
    dF′s = abs(dJ′ - dI):2:(dJ′ + dI)
    I = Int[]
    J = Int[]
    V = Float64[]
    dim = length(basis)
    for q in -1:1
        for (i, dF) in enumerate(dFs)
            idx_offset = offsets[idxFl + i - 1]
            for (i′, dF′) in enumerate(dF′s)
                idx_offset′ = offsets[idxFl′ + i′ - 1]
                for (di, dmF) in enumerate(-dF:2:dF), (di′, dmF′) in enumerate(-dF′:2:dF′)
                    m = dipole_branch(q, dJ, dJ′, dI, dF, dF′, dmF, dmF′)
                    if m == 0
                        continue
                    end
                    push!(I, di + idx_offset)
                    push!(J, di′ + idx_offset′)
                    push!(V, sqrtΓ * m)
                end
            end
            push!(Js, Operator(basis, sparse(I, J, V, dim, dim)))
            empty!(I)
            empty!(J)
            empty!(V)
        end
    end
    return
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
                             sparse([dipole_branch(q, dJ, dJ′, dI, dF, dF′, dmF, dmF′)
                                     for dmF in -dF:2:dF, dmF′ in -dF′:2:dF′]))
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
    op_dagger::O
    nh0::O
    nh0_dagger::O

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
        B = typeof(basis)

        offsets = cumsum(length(b) for b in basis.bases)
        insert!(offsets, 1, 0)

        J = Operator{B,B,SparseMatrixCSC{Float64,Int}}[]
        add_J!(J, basis, offsets, ΓP_S, 1, 1, 1, 1, 3)
        add_J!(J, basis, offsets, ΓP_D, 3, 1, 1, 5, 3)
        add_J!(J, basis, offsets, ΓB_S, 1, 1, 1, 1, 7)
        add_J!(J, basis, offsets, ΓB_D, 3, 1, 1, 5, 7)
        Jdagger = dagger.(J)

        nh0 = h0
        nh0_dagger = copy(h0)

        for (j, jd) in zip(J, Jdagger)
            nhj = jd * j
            nh0 .-= (0.5im) .* nhj
            nh0_dagger .+= (0.5im) .* nhj
        end

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

        return new{B,typeof(h0)}(
            basis, zero(h0), zero(h0), nh0, nh0_dagger, J, Jdagger,
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
    update!(sys, 0.0)
    cb = let update! = update!, sys = sys
        function (t, rho)
            update!(sys, t)
            return (sys.op, sys.op_dagger, sys.J, sys.Jdagger)
        end
    end
    return timeevolution.master_nh_dynamic(range(0, tlen, npoints), ρ0, cb; kws...)
end

function no_drive!(sys::Yb171Sys, t)
    if t == 0
        sys.op .= sys.nh0
        sys.op_dagger .= sys.nh0_dagger
    end
end

struct SingleDrive
    freq::Float64
    pol::NTuple{3,ComplexF64}
end

struct Drives
    d370::Vector{SingleDrive}
    d935::Vector{SingleDrive}
    buff::SparseMatrixCSC{ComplexF64,Int}
    Drives(d370=SingleDrive[], d935=SingleDrive[]) =
        new(d370, d935, spzeros(ComplexF64, 20, 20))
end

function copy_scaled!(tgt::SparseMatrixCSC, src::SparseMatrixCSC, scale)
    resize!(tgt.colptr, length(src.colptr))
    tgt.colptr .= src.colptr
    resize!(tgt.rowval, length(src.rowval))
    tgt.rowval .= src.rowval
    resize!(tgt.nzval, length(src.nzval))
    tgt.nzval .= src.nzval .* scale
    return tgt
end

function add_drive!(sys::Yb171Sys, drives::Drives, dOP, scale)
    if scale == 0
        return
    end
    copy_scaled!(drives.buff, dOP.data, scale)
    sys.op.data .+= drives.buff
    sys.op_dagger.data .+= drives.buff
    return
end

function (drives::Drives)(sys::Yb171Sys, t)
    sys.op .= sys.nh0
    sys.op_dagger .= sys.nh0_dagger

    Ω370 = (zero(ComplexF64), zero(ComplexF64), zero(ComplexF64))
    for d370 in drives.d370
        Ω370 = Ω370 .+ d370.pol .* cis(d370.freq * t)
    end
    add_drive!(sys, drives, sys.dP_Sσ⁻⁺, Ω370[1])
    add_drive!(sys, drives, sys.dP_Sσ⁻⁻, Ω370[1]')
    add_drive!(sys, drives, sys.dP_Sπ⁺, Ω370[1])
    add_drive!(sys, drives, sys.dP_Sπ⁻, Ω370[1]')
    add_drive!(sys, drives, sys.dP_Sσ⁺⁺, Ω370[3])
    add_drive!(sys, drives, sys.dP_Sσ⁺⁻, Ω370[3]')

    Ω935 = (zero(ComplexF64), zero(ComplexF64), zero(ComplexF64))
    for d935 in drives.d935
        Ω935 = Ω935 .+ d935.pol .* cis(d935.freq * t)
    end
    add_drive!(sys, drives, sys.dB_Dσ⁻⁺, Ω935[1])
    add_drive!(sys, drives, sys.dB_Dσ⁻⁻, Ω935[1]')
    add_drive!(sys, drives, sys.dB_Dπ⁺, Ω935[1])
    add_drive!(sys, drives, sys.dB_Dπ⁻, Ω935[1]')
    add_drive!(sys, drives, sys.dB_Dσ⁺⁺, Ω935[3])
    add_drive!(sys, drives, sys.dB_Dσ⁺⁻, Ω935[3]')
    return
end

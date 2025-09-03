#!/usr/bin/julia

using QuantumOptics
using LinearAlgebra
using WignerSymbols
using SparseArrays

using NaCsSim.Master
using NaCsCalc.Atomic: g_sum, g_s, g_l

const g_S1_2 = g_s
const g_P1_2 = g_sum(1/2, 1, g_l, 1/2, g_s)
const g_D3_2 = g_sum(3/2, 2, g_l, 1/2, g_s)
const g_F7_2 = g_sum(7/2, 3, g_l, 1/2, g_s)
const g_B1_2 = g_sum(1/2, 3/2, g_sum(3/2, 7/2, g_F7_2, 2, g_l), 1, g_s)

const μB = 2π * 1.3996244917e6 # Hz/G

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
    V = Float64[]

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

        cg1l = clebschgordan(dJ//2, dmj1//2, dI//2, 1//2, Fl, mF)
        cg2l = clebschgordan(dJ//2, dmj2//2, dI//2, -1//2, Fl, mF)
        cg1h = clebschgordan(dJ//2, dmj1//2, dI//2, 1//2, Fh, mF)
        cg2h = clebschgordan(dJ//2, dmj2//2, dI//2, -1//2, Fh, mF)

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
        wigner6j(dJ//2, dJ′//2, 1, dF′//2, dF//2, dI//2) *
        wigner3j(dF′//2, 1, dF//2, dmF′//2, q, -dmF//2))
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

struct Yb171Data{Basis,O,CO}
    basis::Basis
    nh0_dagger::CO

    J::Vector{O}

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
    function Yb171Data(B)
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

        nh0_dagger = Operator(basis, complex(h0.data))

        for (j, jd) in zip(J, Jdagger)
            nhj = jd * j
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

        return new{B,typeof(h0),typeof(nh0_dagger)}(
            basis, nh0_dagger, J,
            dP_Sσ⁺⁺, dP_Sσ⁺⁻, dP_Sσ⁻⁺, dP_Sσ⁻⁻, dP_Sπ⁺, dP_Sπ⁻,
            dP_Dσ⁺⁺, dP_Dσ⁺⁻, dP_Dσ⁻⁺, dP_Dσ⁻⁻, dP_Dπ⁺, dP_Dπ⁻,
            dB_Sσ⁺⁺, dB_Sσ⁺⁻, dB_Sσ⁻⁺, dB_Sσ⁻⁻, dB_Sπ⁺, dB_Sπ⁻,
            dB_Dσ⁺⁺, dB_Dσ⁺⁻, dB_Dσ⁻⁺, dB_Dσ⁻⁻, dB_Dπ⁺, dB_Dπ⁻
        )
    end
end

const Yb171Sys{Basis,O,CO} = Master.System{Float64,20,Yb171Data{Basis,O,CO}}

function Yb171Sys(B)
    data = Yb171Data(B)
    return Master.System{Float64,20}([j.data for j in data.J], data)
end

function get_ψ(sys::Yb171Sys, name::Symbol, F, mF)
    ψ = Ket(sys.data.basis)
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

function evolve(drive, sys::Yb171Sys, ρ0, tspan; kws...)
    fout(t, x) = Operator(ρ0.basis_l, ρ0.basis_r, copy(x))
    return Master.evolve(drive, sys, ρ0.data, tspan; fout=fout, kws...)
end

function Master.init!(::Nothing, sys::Yb171Sys)
    sys.op_dagger .= sys.data.nh0_dagger.data
end
function Master.update!(::Nothing, sys::Yb171Sys, t)
end

struct SingleDrive
    freq::Float64
    pol::NTuple{3,ComplexF64}
end

struct Drives
    d370::Vector{SingleDrive}
    d935::Vector{SingleDrive}

    nh0_dagger::Vector{ComplexF64}

    dP_Sσ⁺x::Vector{Float64}
    dP_Sσ⁺y::Vector{Float64}
    dP_Sσ⁻x::Vector{Float64}
    dP_Sσ⁻y::Vector{Float64}
    dP_Sπx::Vector{Float64}
    dP_Sπy::Vector{Float64}

    dB_Dσ⁺x::Vector{Float64}
    dB_Dσ⁺y::Vector{Float64}
    dB_Dσ⁻x::Vector{Float64}
    dB_Dσ⁻y::Vector{Float64}
    dB_Dπx::Vector{Float64}
    dB_Dπy::Vector{Float64}

    Drives(d370=SingleDrive[], d935=SingleDrive[]) =
        new(d370, d935, ComplexF64[],
            Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
            Float64[], Float64[], Float64[], Float64[], Float64[], Float64[])
end

function set_vector!(tgt, src)
    resize!(tgt, length(src))
    tgt .= src
    return
end

@inline function add_drive!(op_dagger_nzval, dop_nzvalx, dop_nzvaly, Ω)
    if Ω == 0
        return
    end
    @inbounds @simd ivdep for i in 1:length(op_dagger_nzval)
        v = op_dagger_nzval[i]
        op_dagger_nzval[i] = complex(muladd(dop_nzvalx[i], real(Ω), real(v)),
                                     muladd(dop_nzvaly[i], imag(Ω), imag(v)))
    end
end

function Master.init!(drives::Drives, sys::Yb171Sys)
    data = sys.data
    idx = Ref(1)
    idx_map = Dict{Pair{Int,Int},Int}()
    function add_op(m)
        Master.foreach_nz(m) do row, col, _
            if !haskey(idx_map, row=>col)
                idx_map[row=>col] = idx[]
                idx[] += 1
            end
        end
    end
    add_op(data.nh0_dagger.data)
    has370 = (false, false, false)
    for d370 in drives.d370
        has370 = has370 .| (d370.pol .!= 0)
    end
    if has370[1]
        add_op(data.dP_Sσ⁻⁺.data)
        add_op(data.dP_Sσ⁻⁻.data)
    end
    if has370[2]
        add_op(data.dP_Sπ⁺.data)
        add_op(data.dP_Sπ⁻.data)
    end
    if has370[3]
        add_op(data.dP_Sσ⁺⁺.data)
        add_op(data.dP_Sσ⁺⁻.data)
    end

    has935 = (false, false, false)
    for d935 in drives.d935
        has935 = has935 .| (d935.pol .!= 0)
    end
    if has935[1]
        add_op(data.dB_Dσ⁻⁺.data)
        add_op(data.dB_Dσ⁻⁻.data)
    end
    if has935[2]
        add_op(data.dB_Dπ⁺.data)
        add_op(data.dB_Dπ⁻.data)
    end
    if has935[3]
        add_op(data.dB_Dσ⁺⁺.data)
        add_op(data.dB_Dσ⁺⁻.data)
    end

    I = Int[]
    J = Int[]
    V = Int[]
    for ((i, j), v) in idx_map
        push!(I, i)
        push!(J, j)
        push!(V, v)
    end

    dummy = sparse(I, J, V, 20, 20)
    # Since V is an array of 1:n, sortperm inverse it so that V[order[i]] == i
    order = sortperm(V)

    set_vector!(sys.op_dagger.colptr, dummy.colptr)
    set_vector!(sys.op_dagger.rowval, dummy.rowval)
    resize!(sys.op_dagger.nzval, length(dummy.nzval))

    function remap_values!(m, op, cond=true)
        v = op.data
        resize!(m, length(dummy.nzval))
        if cond
            for i in 1:length(dummy.nzval)
                nzidx = order[dummy.nzval[i]]
                m[i] = v[I[nzidx], J[nzidx]]
            end
        else
            m .= 0
        end
        return
    end

    remap_values!(drives.nh0_dagger, data.nh0_dagger)

    function remap_offdiag_xy!(σx, σy, op⁺, op⁻, cond)
        remap_values!(σx, op⁺, cond)
        remap_values!(σy, op⁻, cond)
        for i in 1:length(σx)
            v⁺ = σx[i]
            v⁻ = σy[i]

            σx[i] = v⁺ + v⁻
            σy[i] = v⁺ - v⁻
        end
    end

    remap_offdiag_xy!(drives.dP_Sσ⁻x, drives.dP_Sσ⁻y,
                      data.dP_Sσ⁻⁺, data.dP_Sσ⁻⁻, has370[1])
    remap_offdiag_xy!(drives.dP_Sπx, drives.dP_Sπy,
                      data.dP_Sπ⁺, data.dP_Sπ⁻, has370[2])
    remap_offdiag_xy!(drives.dP_Sσ⁺x, drives.dP_Sσ⁺y,
                      data.dP_Sσ⁺⁺, data.dP_Sσ⁺⁻, has370[3])

    remap_offdiag_xy!(drives.dB_Dσ⁻x, drives.dB_Dσ⁻y,
                      data.dB_Dσ⁻⁺, data.dB_Dσ⁻⁻, has935[1])
    remap_offdiag_xy!(drives.dB_Dπx, drives.dB_Dπy,
                      data.dB_Dπ⁺, data.dB_Dπ⁻, has935[2])
    remap_offdiag_xy!(drives.dB_Dσ⁺x, drives.dB_Dσ⁺y,
                      data.dB_Dσ⁺⁺, data.dB_Dσ⁺⁻, has935[3])
end
function Master.update!(drives::Drives, sys::Yb171Sys, t)
    op_dagger_nzval = sys.op_dagger.nzval

    nh0_dagger = drives.nh0_dagger
    @inbounds @simd ivdep for i in 1:length(op_dagger_nzval)
        op_dagger_nzval[i] = nh0_dagger[i]
    end

    Ω370 = (zero(ComplexF64), zero(ComplexF64), zero(ComplexF64))
    for d370 in drives.d370
        Ω370 = muladd.(d370.pol, cis(-d370.freq * t), Ω370)
    end
    add_drive!(op_dagger_nzval, drives.dP_Sσ⁻x, drives.dP_Sσ⁻y, Ω370[1])
    add_drive!(op_dagger_nzval, drives.dP_Sπx, drives.dP_Sπy, Ω370[2])
    add_drive!(op_dagger_nzval, drives.dP_Sσ⁺x, drives.dP_Sσ⁺y, Ω370[3])

    Ω935 = (zero(ComplexF64), zero(ComplexF64), zero(ComplexF64))
    for d935 in drives.d935
        Ω935 = muladd.(d935.pol, cis(-d935.freq * t), Ω935)
    end
    add_drive!(op_dagger_nzval, drives.dB_Dσ⁻x, drives.dB_Dσ⁻y, Ω935[1])
    add_drive!(op_dagger_nzval, drives.dB_Dπx, drives.dB_Dπy, Ω935[2])
    add_drive!(op_dagger_nzval, drives.dB_Dσ⁺x, drives.dB_Dσ⁺y, Ω935[3])
    return
end

function print_prob(ρ)
    get_p(i) = round(real(ρ[i, i]) * 100, digits=2)
    println("S0:  0: $(get_p(1))")
    println("S1: -1: $(get_p(2)),  0: $(get_p(3)),  1: $(get_p(4))")
    println("P0:  0: $(get_p(5))")
    println("P1: -1: $(get_p(6)),  0: $(get_p(7)),  1: $(get_p(8))")
    println("D1: -1: $(get_p(9)),  0: $(get_p(10)),  1: $(get_p(11))")
    println("D2: -2: $(get_p(12)), -1: $(get_p(13)),  0: $(get_p(14)),  1: $(get_p(15)),  2: $(get_p(16))")
    println("B0:  0: $(get_p(17))")
    println("B1: -1: $(get_p(18)),  0: $(get_p(19)),  1: $(get_p(20))")
end

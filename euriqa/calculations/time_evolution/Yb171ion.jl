#!/usr/bin/julia

using QuantumOptics
using LinearAlgebra
using CGcoefficient
using SparseArrays

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

function foreach_nz(cb, m)
    for col in 1:size(m, 2)
        for r in nzrange(m, col)
            row = rowvals(m)[r]
            v = nonzeros(m)[r]
            cb(row, col, v)
        end
    end
end

# State order,
# S1/2 (F=0, F=1), P1/2 (F=0, F=1), D3/2 (F=1, F=2), [3/2]1/2 (F=0, F=1)
#         1,   2,          3,   4,          5,   6,              7,   8

struct LinearOP{Map}
end

@inline function map_op!(::LinearOP{Map}, tgt, src) where Map
    @inbounds for (yi, yo, xi, xo, v) in Map
        tgt[xo, yo] += src[xi, yi] * v
    end
end

struct Yb171Sys{Basis,O,CO,JMap}
    basis::Basis
    op_dagger::CO
    nh0::CO
    nh0_dagger::CO

    J::Vector{O}
    Jdagger::Vector{O}
    Jop::LinearOP{JMap}

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

        jmap = Dict{NTuple{4,Int},Float64}()
        function add_jterm(xi, yi, xo, yo, v)
            key = yi, yo, xi, xo
            if haskey(jmap, key)
                jmap[key] += v
            else
                jmap[key] = v
            end
            return
        end
        for j in J
            foreach_nz(j.data) do x1, y1, v1
                foreach_nz(j.data) do x2, y2, v2
                    add_jterm(y1, y2, x1, x2, v1 * v2)
                end
            end
        end
        jop = Tuple{Int,Int,Int,Int,Float64}[]
        for ((yi, yo, xi, xo), v) in jmap
            push!(jop, (yi, yo, xi, xo, v))
        end
        sort!(jop)
        jop = tuple(jop...)

        nh0 = Operator(basis, complex(h0.data))
        nh0_dagger = copy(nh0)

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

        return new{B,typeof(h0),typeof(nh0),jop}(
            basis, zero(nh0), nh0, nh0_dagger, J, Jdagger, LinearOP{jop}(),
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

function init! end
function update! end

@inline function do_mul!(result, B, M::SparseMatrixCSC)
    nrow = size(result, 1)
    @inbounds prev_colptr = M.colptr[1]
    @inbounds for col in 1:M.n
        filled = false
        colptr = M.colptr[col + 1]
        for i in prev_colptr:colptr - 1
            val = M.nzval[i]
            row = M.rowval[i]
            if filled
                @simd ivdep for j in 1:nrow
                    result[j, col] += val * B[j, row]
                end
            else
                @simd ivdep for j in 1:nrow
                    result[j, col] = val * B[j, row]
                end
            end
            filled = true
        end
        prev_colptr = colptr
        if !filled
            @simd ivdep for j in 1:nrow
                result[j, col] = 0
            end
        end
    end
end

function evolve(drive, sys::Yb171Sys, ρ0, tlen, npoints=1001; kws...)
    init!(drive, sys)
    function dmaster_(t, rho, drho)
        rho_data = rho.data
        drho_data = drho.data

        update!(drive, sys, t)
        nrow = size(drho_data, 1)
        do_mul!(drho_data, rho_data, sys.op_dagger.data)

        # compute -i * drho^dagger + i * drho
        @inbounds for i in 1:nrow
            drho_data[i, i] = -2 * imag(drho_data[i, i])
            for j in (i + 1):nrow
                v1 = drho_data[j, i]
                v2 = drho_data[i, j]

                re_o = -imag(v1) - imag(v2)
                im_o = real(v2) - real(v1)

                drho_data[i, j] = complex(re_o, im_o)
                drho_data[j, i] = complex(re_o, -im_o)
            end
        end

        map_op!(sys.Jop, drho_data, rho_data)
        return drho
    end
    return timeevolution.integrate_master(range(0, tlen, npoints), dmaster_,
                                          ρ0, nothing; kws...)
end

function init!(::Nothing, sys::Yb171Sys)
    sys.op_dagger .= sys.nh0_dagger
end
function update!(::Nothing, sys::Yb171Sys, t)
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
        op_dagger_nzval[i] += complex(dop_nzvalx[i] * real(Ω),
                                      dop_nzvaly[i] * imag(Ω))
    end
end

function init!(drives::Drives, sys::Yb171Sys)
    idx = Ref(1)
    idx_map = Dict{Pair{Int,Int},Int}()
    function add_op(m)
        foreach_nz(m) do row, col, _
            if !haskey(idx_map, row=>col)
                idx_map[row=>col] = idx[]
                idx[] += 1
            end
        end
    end
    add_op(sys.nh0_dagger.data)
    has370 = (false, false, false)
    for d370 in drives.d370
        has370 = has370 .| (d370.pol .!= 0)
    end
    if has370[1]
        add_op(sys.dP_Sσ⁻⁺.data)
        add_op(sys.dP_Sσ⁻⁻.data)
    end
    if has370[2]
        add_op(sys.dP_Sπ⁺.data)
        add_op(sys.dP_Sπ⁻.data)
    end
    if has370[3]
        add_op(sys.dP_Sσ⁺⁺.data)
        add_op(sys.dP_Sσ⁺⁻.data)
    end

    has935 = (false, false, false)
    for d935 in drives.d935
        has935 = has935 .| (d935.pol .!= 0)
    end
    if has935[1]
        add_op(sys.dB_Dσ⁻⁺.data)
        add_op(sys.dB_Dσ⁻⁻.data)
    end
    if has935[2]
        add_op(sys.dB_Dπ⁺.data)
        add_op(sys.dB_Dπ⁻.data)
    end
    if has935[3]
        add_op(sys.dB_Dσ⁺⁺.data)
        add_op(sys.dB_Dσ⁺⁻.data)
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

    set_vector!(sys.op_dagger.data.colptr, dummy.colptr)
    set_vector!(sys.op_dagger.data.rowval, dummy.rowval)
    resize!(sys.op_dagger.data.nzval, length(dummy.nzval))

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

    remap_values!(drives.nh0_dagger, sys.nh0_dagger)

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
                      sys.dP_Sσ⁻⁺, sys.dP_Sσ⁻⁻, has370[1])
    remap_offdiag_xy!(drives.dP_Sπx, drives.dP_Sπy,
                      sys.dP_Sπ⁺, sys.dP_Sπ⁻, has370[2])
    remap_offdiag_xy!(drives.dP_Sσ⁺x, drives.dP_Sσ⁺y,
                      sys.dP_Sσ⁺⁺, sys.dP_Sσ⁺⁻, has370[3])

    remap_offdiag_xy!(drives.dB_Dσ⁻x, drives.dB_Dσ⁻y,
                      sys.dB_Dσ⁻⁺, sys.dB_Dσ⁻⁻, has935[1])
    remap_offdiag_xy!(drives.dB_Dπx, drives.dB_Dπy,
                      sys.dB_Dπ⁺, sys.dB_Dπ⁻, has935[2])
    remap_offdiag_xy!(drives.dB_Dσ⁺x, drives.dB_Dσ⁺y,
                      sys.dB_Dσ⁺⁺, sys.dB_Dσ⁺⁻, has935[3])
end
function update!(drives::Drives, sys::Yb171Sys, t)
    op_dagger_nzval = sys.op_dagger.data.nzval

    nh0_dagger = drives.nh0_dagger
    @inbounds @simd ivdep for i in 1:length(op_dagger_nzval)
        op_dagger_nzval[i] = nh0_dagger[i]
    end

    Ω370 = (zero(ComplexF64), zero(ComplexF64), zero(ComplexF64))
    for d370 in drives.d370
        Ω370 = Ω370 .+ d370.pol .* cis(-d370.freq * t)
    end
    add_drive!(op_dagger_nzval, drives.dP_Sσ⁻x, drives.dP_Sσ⁻y, Ω370[1])
    add_drive!(op_dagger_nzval, drives.dP_Sπx, drives.dP_Sπy, Ω370[2])
    add_drive!(op_dagger_nzval, drives.dP_Sσ⁺x, drives.dP_Sσ⁺y, Ω370[3])

    Ω935 = (zero(ComplexF64), zero(ComplexF64), zero(ComplexF64))
    for d935 in drives.d935
        Ω935 = Ω935 .+ d935.pol .* cis(-d935.freq * t)
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

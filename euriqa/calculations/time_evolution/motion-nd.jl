#!/usr/bin/julia

using QuantumOptics
using StaticArrays
using PyPlot
using SparseArrays

using NaCsSim.Master

function motion_energy(ωs, ns)
    return sum(ωs .* ns)
end

function _collect_motion_states!(res, ωs::NTuple{N,T}, lb, ub, prefix) where {N,T}
    idx = length(prefix) + 1
    ω = ωs[idx]
    ubn = floor(Int, ub / ω)
    if idx == N
        for n in max(ceil(Int, lb / ω), 0):ubn
            push!(res, (prefix..., n))
        end
        return
    end
    for n in 0:ubn
        de = n * ω
        _collect_motion_states!(res, ωs, lb - de, ub - de, (prefix..., n))
    end
end

function collect_motion_states(ωs::NTuple{N,T}, lb, ub) where {N,T}
    res = NTuple{N,Int}[]
    _collect_motion_states!(res, ωs, lb, ub, ())
    return res
end

function displace_sparse(states, αs::NTuple{N,T}) where {N,T}
    nrow = length(states)
    I = Int[]
    J = Int[]
    V = T[]
    for i in 1:nrow
        for j in i:nrow
            v = prod(displace_analytical.(αs, states[i], states[j]))
            if abs(v) >= 1e-10
                push!(I, i)
                push!(J, j)
                push!(V, v)
                if i != j
                    push!(I, j)
                    push!(J, i)
                    push!(V, prod(displace_analytical.(αs, states[j], states[i])))
                end
            end
        end
    end
    return sparse(I, J, V, nrow, nrow)
end

struct MotionNDData{N}
    E::Float64
    ωs::NTuple{N,Float64}
    ηs::NTuple{N,Float64}
    states::Vector{NTuple{N,Int}}
    disp::SparseMatrixCSC{Float64,Int}
    disp_dagger::SparseMatrixCSC{Float64,Int}
end
function MotionNDData(E, ωs::NTuple{N}, ηs::NTuple{N}, states) where N
    disp = displace_sparse(states, ηs)
    return MotionNDData{N}(E, ωs, ηs, states, disp, disp')
end
const MotionND{N,nrow} = Master.SystemCoherent{Float64,nrow,MotionNDData{N}}

function Master.init!(drive, sys::MotionND)
    data = sys.data
    disp = data.disp
    nmotion = length(data.states)
    disp_nnz = length(disp.rowval)
    half_op_nnz = disp_nnz + nmotion
    op_nnz = half_op_nnz * 2

    colptr = sys.op.colptr
    resize!(colptr, 2 * nmotion + 1)
    colptr[1] = 1
    rowval = sys.op.rowval
    resize!(rowval, op_nnz)
    nzval = sys.op.nzval
    resize!(nzval, op_nnz)

    prev_disp_colptr = 1
    prev_colptr1 = 1
    for i in 1:nmotion
        disp_colptr = disp.colptr[i + 1]
        colptr1 = disp_colptr + i
        colptr[i + 1] = colptr1
        colptr[i + 1 + nmotion] = colptr1 + half_op_nnz

        E_m = motion_energy(data.ωs, data.states[i])

        # region to fill in rowval:
        # [prev_colptr1, colptr1)
        # [prev_colptr1 + half_op_nnz, colptr1 + half_op_nnz)
        # for [prev_colptr1, colptr1), the first element is the diagonal one
        rowval[prev_colptr1] = i
        nzval[prev_colptr1] = E_m - data.E / 2
        # for [prev_colptr1 + half_op_nnz, colptr1 + half_op_nnz),
        # the last element is the diagonal one
        rowval[colptr1 + half_op_nnz - 1] = i + nmotion
        nzval[colptr1 + half_op_nnz - 1] = E_m + data.E / 2
        for disp_cp in prev_disp_colptr:disp_colptr - 1
            row = disp.rowval[disp_cp]
            rowval[disp_cp + i] = row + nmotion
            rowval[disp_cp + i - 1 + half_op_nnz] = row
        end

        prev_disp_colptr = disp_colptr
        prev_colptr1 = colptr1
    end
end

function Master.update!(drive, sys::MotionND, t)
    Ω = drive(t)

    data = sys.data
    disp = data.disp
    disp_dagger = data.disp_dagger
    nmotion = length(data.states)
    disp_nnz = length(disp.rowval)
    half_op_nnz = disp_nnz + nmotion

    nzval = sys.op.nzval

    prev_disp_colptr = 1
    prev_colptr1 = 1
    for i in 1:nmotion
        disp_colptr = disp.colptr[i + 1]
        colptr1 = disp_colptr + i

        for disp_cp in prev_disp_colptr:disp_colptr - 1
            v = disp.nzval[disp_cp]
            nzval[disp_cp + i] = v * Ω
            v_dagger = disp_dagger.nzval[disp_cp]
            nzval[disp_cp + i - 1 + half_op_nnz] = v_dagger * conj(Ω)
        end
        prev_disp_colptr = disp_colptr
        prev_colptr1 = colptr1
    end
    return
end

struct Drive{T,F}
    Ω_2::T
    f::F
end
(drive::Drive{T,Nothing} where T)(t) = drive.Ω_2
@inline (drive::Drive)(t) = drive.Ω_2 * drive.f(t)

function evolve(E, ωs, ηs, ψs0, n0s, Ω, tlen, npoints=1001;
                δmax=0, Ωprofile::F=nothing, kws...) where F
    E0 = motion_energy(ωs, n0s)
    dE = hypot(Ω, δmax) * 10
    mstates = collect_motion_states(ωs, max(E0 - dE, 0.0), E0 + dE)
    N = length(mstates)
    data = MotionNDData(E, ωs, ηs, mstates)
    sys = Master.SystemCoherent{Float64,nothing}(data, 2 * N)
    m_basis = NLevelBasis(N)
    ψ0 = nlevelstate(m_basis, findfirst(==(n0s), mstates)) ⊗ ψs0
    basis = m_basis ⊗ SpinBasis(1//2)
    ts, xs = Master.evolve(Drive(Ω / 2, Ωprofile), sys, ψ0.data, tlen, npoints; kws...)
    return ts, [Ket(basis, x) for x in xs], sys
end

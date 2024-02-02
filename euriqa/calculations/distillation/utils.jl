#!/usr/bin/julia

using LinearAlgebra

function expand_operator(::Val{N}, O, op_dim) where N
    indices = CartesianIndices(ntuple(_->2, N))
    O_lindices = LinearIndices(ntuple(_->2, length(op_dim)))
    M_sz = length(indices)

    M = Matrix{eltype(O)}(undef, M_sz, M_sz)

    other_dim = [i for i in 1:N if !(i in op_dim)]

    for (li1, ci1) in zip(1:M_sz, indices)
        mi1 = O_lindices[getindex.(ci1, op_dim)...]
        for (li2, ci2) in zip(1:M_sz, indices)
            is_diag = true
            for i in other_dim
                if ci1[i] != ci2[i]
                    is_diag = false
                    break
                end
            end
            if is_diag
                mi2 = O_lindices[getindex.(ci2, op_dim)...]
                M[li2, li1] = O[mi2, mi1]
            else
                M[li2, li1] = 0
            end
        end
    end
    return M
end

function CnZ(n, εz)
    O = Matrix{ComplexF64}(I, 2^n, 2^n)
    O[end, end] = -cis(εz * π)
    return O
end

function CNOT(εr)
    C = -sinpi(εr / 2)
    S = cospi(εr / 2)
    return [1 0 0 0
            0 1 0 0
            0 0 C S
            0 0 S C]
end

function SWAP(εr)
    C = -sinpi(εr / 2)
    S = cospi(εr / 2)
    return [1 0 0 0
            0 C S 0
            0 S C 0
            0 0 0 1]
end

const HGate = [1 1
               1 -1] ./ sqrt(2)
const XGate = [0 1
               1 0]

function generation_circuit(εr, εz)
    N = 4
    VN = Val(N)
    U = expand_operator(VN, HGate, (2,)) * expand_operator(VN, HGate, (4,))
    U = CnZ(4, εz) * U
    U = (expand_operator(VN, HGate, (1,)) *
        expand_operator(VN, HGate * XGate, (3,)) *
        expand_operator(VN, HGate, (4,)) * U)
    U = expand_operator(VN, CnZ(2, εz), (1, 3)) * U
    return expand_operator(VN, HGate, (3,)) * U
end

function distillation_circuit(εr, εz)
    N = 4
    VN = Val(N)
    return (expand_operator(VN, SWAP(εr) * CNOT(εr), (1, 2)) *
        expand_operator(VN, SWAP(εr) * CNOT(εr), (3, 4)))
end

measure_bit(::Val{N}, ψ::AbstractVector, bit) where N =
    measure_bit(Val(N), ψ * ψ', bit)
function measure_bit(::Val{N}, ρ::AbstractMatrix, bit) where N
    ρ0 = similar(ρ)
    ρ1 = similar(ρ)
    indices = CartesianIndices(ntuple(_->2, N))
    M_sz = length(indices)
    for (li1, ci1) in zip(1:M_sz, indices)
        v1 = ci1[bit] == 2
        for (li2, ci2) in zip(1:M_sz, indices)
            v2 = ci2[bit] == 2
            if v1 != v2
                ρ0[li2, li1] = 0
                ρ1[li2, li1] = 0
            elseif v1
                ρ0[li2, li1] = 0
                ρ1[li2, li1] = ρ[li2, li1]
            else
                ρ0[li2, li1] = ρ[li2, li1]
                ρ1[li2, li1] = 0
            end
        end
    end
    p0 = tr(ρ0)
    p1 = tr(ρ1)
    if p0 != 0
        ρ0 ./= p0
    end
    if p1 != 0
        ρ1 ./= p1
    end
    return [(p0, ρ0), (p1, ρ1)]
end

const ψ0 = zeros(16)
ψ0[6] = 1

function fidelity(ψ::AbstractVector, ψ0)
    return abs2(dot(ψ, ψ0))
end

function fidelity(ρ::AbstractMatrix, ψ0)
    return tr(ψ0' * ρ * ψ0)
end

fidelity_per_bit(a, b) = sqrt(fidelity(a, b))

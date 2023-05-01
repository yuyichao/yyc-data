#!/usr/bin/julia

using LinearAlgebra

struct System{T}
    N::Int
    M::Matrix{T}
    System{T}(N) where T = new{T}(N, zeros(T, (N^2, N^2)))
end

function propagate(sys::System{T}, ρ0, t) where T
    ρ = similar(ρ0, T)
    mul!(vec(ρ), exp(sys.M .* t), vec(ρ0))
    return ρ
end

# The system state vector v is the flattened density matrix.
# (1-based indexing)
# v(i + (j - 1) * N) = ρ(i, j)
function add_H!(sys::System{T}, H) where T
    N = sys.N
    M = sys.M
    # d ρᵢⱼ/dt = 1/im * (Hρ - ρH)ᵢⱼ
    #          = 1/im * ∑ₖ(Hᵢₖρₖⱼ-ρᵢₖHₖⱼ)
    for i in 1:N
        for j in 1:N
            mi = i + (j - 1) * N
            for k in 1:N
                M[mi, k + (j - 1) * N] += H[i, k] / im
                M[mi, i + (k - 1) * N] -= H[k, j] / im
            end
        end
    end
    return sys
end
function add_C!(sys::System{T}, C) where T
    N = sys.N
    M = sys.M
    L = C * C'
    # d ρᵢⱼ/dt = (C ρ C†)ᵢⱼ - (C C† ρ + ρ C C†)ᵢⱼ / 2
    #          = ∑ₖₗ Cᵢₖ Cⱼₗ* ρₖₗ  - ∑ₖ(Lᵢₖ ρₖⱼ + ρᵢₖ Lₖⱼ) / 2
    for i in 1:N
        for j in 1:N
            mi = i + (j - 1) * N
            for k in 1:N
                for l in 1:N
                    M[mi, k + (l - 1) * N] += C[i, k] * conj(C[j, l])
                end
                M[mi, k + (j - 1) * N] -= L[i, k] / 2
                M[mi, i + (k - 1) * N] -= L[k, j] / 2
            end
        end
    end
    return sys
end

new_H(sys::System{T}) where T = zeros(T, (sys.N, sys.N))
function H_add_E!(H, E, i)
    H[i, i] += E
    return H
end

function H_add_Ω!(H, Ω, i, j)
    H[i, j] += Ω
    H[j, i] += conj(Ω)
    return H
end

new_C(sys::System{T}) where T = zeros(T, (sys.N, sys.N))

function C_add_decay!(C, Γ, to, from)
    C[to, from] = sqrt(Γ)
    return C
end

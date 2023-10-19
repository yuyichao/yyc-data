#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Trap
using NaCsPlot
using PyPlot
using LinearAlgebra

function get_H(N, ω_m, δ, Ω, η)
    H = zeros(ComplexF64, 2N, 2N)
    for i in 1:N
        n = i - 1
        H[i, i] += n * ω_m # ground
        H[i + N, i + N] += n * ω_m - δ # excited
    end
    for i in 1:N
        ng = i - 1 # ground
        for j in 1:N
            ne = j - 1 # excited
            M = Ω * Trap.sideband(ng, ne, η) / 2 * (im^mod(abs(j - i), 4))
            H[i, j + N] = M
            H[j + N, i] = conj(M)
        end
    end
    return H
end

function get_ψ(N, n, e=false)
    @assert n < N
    ψ = zeros(ComplexF64, 2N)
    i = n + 1 + (e ? N : 0)
    ψ[i] = 1
    return ψ
end

function get_pe(ψ)
    p = 0.0
    N = length(ψ) ÷ 2
    for i in 1:N
        p += abs2(ψ[i + N])
    end
    return p
end

function get_pg(ψ)
    p = 0.0
    N = length(ψ) ÷ 2
    for i in 1:N
        p += abs2(ψ[i])
    end
    return p
end

mutable struct ExpCache{T,TH,TU}
    const H::TH
    const U::TU
    t::T
    function ExpCache(H::AbstractMatrix{T}) where T
        TR = real(T)
        TC = complex(TR)
        TH = typeof(H)
        TU = Matrix{TC}
        return new{TR,TH,TU}(H, TU(I, size(H)...), 0)
    end
end

function get_U(cache::ExpCache, t)
    if cache.t != t
        cache.t = t
        cache.U .= exp((im * t) .* cache.H)
    end
    return cache.U
end

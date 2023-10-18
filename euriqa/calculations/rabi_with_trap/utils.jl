#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Trap
using NaCsPlot
using PyPlot

function get_H(N, ω_m, δ, Ω, η)
    H = zeros(ComplexF64, 2N, 2N)
    for i in 1:N
        n = i - 1
        H[i, i] += n * ω # ground
        H[i + N, i + N] += n * ω - δ # excited
    end
    for i in 1:N
        ng = i - 1 # ground
        for j in 1:N
            ne = j - 1 # excited
            M = Ω * Trap.sideband(ng, ne, η) / 2
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

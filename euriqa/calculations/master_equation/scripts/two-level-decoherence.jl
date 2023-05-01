#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

const sys = System{ComplexF64}(2)
const H = new_H(sys)
H_add_E!(H, 2π * 40, 2)
H_add_Ω!(H, 2π * 6, 1, 2)
add_H!(sys, H)

const C1 = new_C(sys)
C_add_decay!(C1, 2π * 8.0, 1, 1)
add_C!(sys, C1)

const C2 = new_C(sys)
C_add_decay!(C2, 2π * 8.0, 2, 2)
add_C!(sys, C2)

const ts = range(0, 2, 10001)
const ρ0 = [1 0; 0 0]

function analytical_nodecay(t)
    Ω = 2π * 6
    Δ = 2π * 40
    Ω² = Ω^2
    Ωg² = Ω² + Δ^2
    Ωg = sqrt(Ωg²)
    return Ω² / Ωg² * sin(Ωg * t / 2)^2
end

function approx_decay(t)
    Ω = 2π * 6
    Δ = 2π * 40
    Γ = 2π * 8
    Ω² = Ω^2
    Ωg² = Ω² + Δ^2

    γ = Γ * Ω² / Ωg²
    return 0.5 * (1 - exp(-t * γ))
end

const prefix = joinpath(@__DIR__, "../imgs/two-level-decoherence")

figure()
plot(ts, (t->(real(propagate(sys, ρ0, t)[2, 2]))).(ts))
plot(ts, analytical_nodecay.(ts))
plot(ts, approx_decay.(ts), "--")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

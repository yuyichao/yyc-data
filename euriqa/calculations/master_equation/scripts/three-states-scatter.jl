#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

function get_sys_excited(Ω, Δ, Γ)
    sys = System{ComplexF64}(3)

    H = new_H(sys)
    H_add_E!(H, Δ, 2)
    H_add_Ω!(H, Ω, 1, 2)
    add_H!(sys, H)

    C1 = new_C(sys)
    C_add_decay!(C1, Γ / 2, 1, 2)
    add_C!(sys, C1)

    C2 = new_C(sys)
    C_add_decay!(C2, Γ / 2, 3, 2)
    add_C!(sys, C2)
end

function get_sys_ground(Ω, Δ, Γ)
    sys = System{ComplexF64}(3)

    γ = Γ * Ω^2 / (Γ^2 + 4 * Δ^2 + 2 * Ω^2)

    C1 = new_C(sys)
    C_add_decay!(C1, γ / 2, 1, 1)
    add_C!(sys, C1)

    C2 = new_C(sys)
    C_add_decay!(C2, γ / 2, 3, 1)
    add_C!(sys, C2)
end

const sys_g0 = get_sys_ground(2π * 10, 2π * 1000, 2π * 40)
const sys_e0 = get_sys_excited(2π * 10, 2π * 1000, 2π * 40)
const sys_g1 = get_sys_ground(2π * 10, 2π * 300, 2π * 40)
const sys_e1 = get_sys_excited(2π * 10, 2π * 300, 2π * 40)

const ts = range(0, 400, 10001)
const ρ0 = [1 0 0; 0 0 0; 0 0 0]

const prefix = joinpath(@__DIR__, "../imgs/three-states-scatter")

figure()
plot(ts, (t->(real(propagate(sys_e0, ρ0, t)[1, 1]))).(ts), "C0")
plot(ts, (t->(real(propagate(sys_g0, ρ0, t)[1, 1]))).(ts), "C2--")
plot(ts, (t->(real(propagate(sys_e1, ρ0, t)[1, 1]))).(ts), "C1")
plot(ts, (t->(real(propagate(sys_g1, ρ0, t)[1, 1]))).(ts), "C3--")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

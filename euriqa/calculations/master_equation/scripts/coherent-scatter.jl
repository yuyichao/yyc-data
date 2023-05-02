#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

function get_sys_independent(Ω, Δ, γ)
    sys = System{ComplexF64}(2)

    H = new_H(sys)
    H_add_E!(H, Δ / 2, 1)
    H_add_E!(H, -Δ / 2, 2)
    H_add_Ω!(H, Ω, 1, 2)
    add_H!(sys, H)

    C = new_C(sys)
    C_add_decay!(C, γ / 2, 1, 1)
    add_C!(sys, C)

    C = new_C(sys)
    C_add_decay!(C, γ / 2, 1, 2)
    add_C!(sys, C)

    C = new_C(sys)
    C_add_decay!(C, γ / 2, 2, 1)
    add_C!(sys, C)

    C = new_C(sys)
    C_add_decay!(C, γ / 2, 2, 2)
    add_C!(sys, C)

    return sys
end

function get_sys_coherent(Ω, Δ, γ)
    sys = System{ComplexF64}(2)

    H = new_H(sys)
    H_add_E!(H, Δ / 2, 1)
    H_add_E!(H, -Δ / 2, 2)
    H_add_Ω!(H, Ω, 1, 2)
    add_H!(sys, H)

    C = new_C(sys)
    C_add_decay!(C, γ / 2, 1, 1)
    C_add_decay!(C, γ / 2, 1, 2)
    add_C!(sys, C)

    C = new_C(sys)
    C_add_decay!(C, γ / 2, 2, 1)
    C_add_decay!(C, γ / 2, 2, 2)
    add_C!(sys, C)

    return sys
end

const sys_g0 = get_sys_independent(2π * 0.1, -2π * 0.2, 2π * 0.02)
const sys_e0 = get_sys_coherent(2π * 0.1, -2π * 0.2, 2π * 0.02)
const sys_g1 = get_sys_independent(2π * 0.1, 2π * 0.2, 2π * 0.02)
const sys_e1 = get_sys_coherent(2π * 0.1, 2π * 0.2, 2π * 0.02)

const ts = range(0, 40, 10001)
const dets = 2π * range(-0.5, 0.5, 10001)
const ρ0 = [1 0; 0 0]

function get_22_independent(Ω, Δ, γ, t)
    sys = get_sys_independent(Ω, Δ, γ)
    return real(propagate(sys, ρ0, t)[2, 2])
end

function get_22_coherent(Ω, Δ, γ, t)
    sys = get_sys_coherent(Ω, Δ, γ)
    return real(propagate(sys, ρ0, t)[2, 2])
end

const prefix = joinpath(@__DIR__, "../imgs/coherent-scatter")

figure()
plot(ts, (t->(real(propagate(sys_e0, ρ0, t)[2, 2]))).(ts), "C0")
plot(ts, (t->(real(propagate(sys_g0, ρ0, t)[2, 2]))).(ts), "C2--")
plot(ts, (t->(real(propagate(sys_e1, ρ0, t)[2, 2]))).(ts), "C1")
plot(ts, (t->(real(propagate(sys_g1, ρ0, t)[2, 2]))).(ts), "C3--")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save(prefix)

figure()
plot(dets, get_22_coherent.(2π * 0.1, dets, 2π * 0.02, 5), "C0")
plot(dets, get_22_coherent.(2π * 0.1, dets, 2π * 0.02, 15), "C1")
plot(dets, get_22_coherent.(2π * 0.1, dets, 2π * 0.02, 35), "C2")
plot(dets, get_22_coherent.(2π * 0.1, dets, 2π * 0.02, 75), "C3")
plot(dets, get_22_independent.(2π * 0.1, dets, 2π * 0.02, 5), "C0--")
plot(dets, get_22_independent.(2π * 0.1, dets, 2π * 0.02, 15), "C1--")
plot(dets, get_22_independent.(2π * 0.1, dets, 2π * 0.02, 35), "C2--")
plot(dets, get_22_independent.(2π * 0.1, dets, 2π * 0.02, 75), "C3--")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_det")

NaCsPlot.maybe_show()

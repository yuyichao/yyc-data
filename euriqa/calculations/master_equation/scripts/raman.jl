#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

function get_sys(det)
    sys = System{ComplexF64}(3)
    H = new_H(sys)
    H_add_E!(H, 2π * 1000, 2)
    H_add_E!(H, det, 3)
    H_add_Ω!(H, 2π * 10, 1, 2)
    H_add_Ω!(H, 2π * 10, 2, 3)
    add_H!(sys, H)
    return sys
end

const sys0 = get_sys(0)

const ts = range(0, 40, 10001)
const dets = 2π * range(-0.25, 0.25, 10001)
const ρ0 = [1 0 0; 0 0 0; 0 0 0]

function get_33(det, t)
    sys = get_sys(det)
    return real(propagate(sys, ρ0, t)[3, 3])
end

function analytical(det, t)
    Ω₁ = 2π * 10
    Ω₂ = 2π * 10
    Δ = 2π * 1000
    ΩR = Ω₁ * Ω₂ / 2 / Δ
    ΩR² = ΩR^2
    ΩRg² = ΩR² + det^2
    ΩRg = sqrt(ΩRg²)
    return ΩR² / ΩRg² * sin(ΩRg * t / 2)^2
end

const prefix = joinpath(@__DIR__, "../imgs/raman")

figure()
plot(ts, (t->(real(propagate(sys0, ρ0, t)[1, 1]))).(ts), "C0")
plot(ts, (t->(real(propagate(sys0, ρ0, t)[2, 2]))).(ts), "C4")
plot(ts, (t->(real(propagate(sys0, ρ0, t)[3, 3]))).(ts), "C1")
plot(ts, 1 .- analytical.(0, ts), "C2--")
plot(ts, analytical.(0, ts), "C3--")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save(prefix)

figure()
plot(dets ./ (2π), get_33.(dets, 10), "C0")
plot(dets ./ (2π), analytical.(dets, 10), "C2--")
ylim([0, 1])
grid()
xlabel("det")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_det")

NaCsPlot.maybe_show()

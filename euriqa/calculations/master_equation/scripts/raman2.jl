#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

function get_sys(det, Ω₁, Ω₂)
    sys = System{ComplexF64}(3)
    H = new_H(sys)
    H_add_E!(H, 2π * 1000, 2)
    H_add_E!(H, det, 3)
    H_add_Ω!(H, Ω₁, 1, 2)
    H_add_Ω!(H, Ω₂, 2, 3)
    add_H!(sys, H)
    return sys
end

const dets = 2π * range(-0.8, 0.8, 10001)
const ρ0 = [1 0 0; 0 0 0; 0 0 0]

function get_33(det, Ω₁, Ω₂, t)
    sys = get_sys(det, Ω₁, Ω₂)
    return real(propagate(sys, ρ0, t)[3, 3])
end

function analytical(det, Ω₁, Ω₂, t)
    Δ = 2π * 1000
    ΩR = Ω₁ * Ω₂ / 2 / Δ

    det = det + (Ω₁^2 - Ω₂^2) / 4 / Δ

    ΩR² = ΩR^2
    ΩRg² = ΩR² + det^2
    ΩRg = sqrt(ΩRg²)
    return ΩR² / ΩRg² * sin(ΩRg * t / 2)^2
end

const prefix = joinpath(@__DIR__, "../imgs/raman2")

figure()
plot(dets ./ (2π), get_33.(dets, 2π * 10 * 0.2, 2π * 10 / 0.2, 10), "C0")
plot(dets ./ (2π), get_33.(dets, 2π * 10 * 0.3, 2π * 10 / 0.3, 10), "C1")
plot(dets ./ (2π), get_33.(dets, 2π * 10 * 0.5, 2π * 10 / 0.5, 10), "C2")
plot(dets ./ (2π), get_33.(dets, 2π * 10, 2π * 10, 10), "C3")
plot(dets ./ (2π), get_33.(dets, 2π * 10 / 0.5, 2π * 10 * 0.5, 10), "C4")
plot(dets ./ (2π), get_33.(dets, 2π * 10 / 0.3, 2π * 10 * 0.3, 10), "C5")
plot(dets ./ (2π), get_33.(dets, 2π * 10 / 0.2, 2π * 10 * 0.2, 10), "C6")

plot(dets ./ (2π), analytical.(dets, 2π * 10 * 0.2, 2π * 10 / 0.2, 10), "C2:")
plot(dets ./ (2π), analytical.(dets, 2π * 10 * 0.3, 2π * 10 / 0.3, 10), "C3:")
plot(dets ./ (2π), analytical.(dets, 2π * 10 * 0.5, 2π * 10 / 0.5, 10), "C4:")
plot(dets ./ (2π), analytical.(dets, 2π * 10, 2π * 10, 10), "C5:")
plot(dets ./ (2π), analytical.(dets, 2π * 10 / 0.5, 2π * 10 * 0.5, 10), "C6:")
plot(dets ./ (2π), analytical.(dets, 2π * 10 / 0.3, 2π * 10 * 0.3, 10), "C7:")
plot(dets ./ (2π), analytical.(dets, 2π * 10 / 0.2, 2π * 10 * 0.2, 10), "C8:")

ylim([0, 1])
grid()
xlabel("det")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_det")

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

function get_sys(det)
    sys = System{ComplexF64}(3)

    Δ = 2π * 1000
    Ω = 2π * 10
    Γ = 2π * 40

    H = new_H(sys)
    H_add_E!(H, Δ, 2)
    H_add_E!(H, det / 2, 3)
    H_add_E!(H, -det / 2, 1)
    H_add_Ω!(H, Ω, 1, 2)
    H_add_Ω!(H, Ω, 2, 3)
    add_H!(sys, H)

    C1 = new_C(sys)
    C_add_decay!(C1, Γ / 2, 1, 2)
    add_C!(sys, C1)

    C2 = new_C(sys)
    C_add_decay!(C2, Γ / 2, 3, 2)
    add_C!(sys, C2)

    # γ = Γ * Ω^2 / (Γ^2 + 4 * Δ^2 + 2 * Ω^2)

    # C = new_C(sys)
    # C_add_decay!(C, γ / 2, 1, 1)
    # C_add_decay!(C, γ / 2, 1, 3)
    # add_C!(sys, C)

    # C = new_C(sys)
    # C_add_decay!(C, γ / 2, 3, 1)
    # C_add_decay!(C, γ / 2, 3, 3)
    # add_C!(sys, C)

    return sys
end

const sys0 = get_sys(0)

const ts = range(0, 400, 10001)
const dets = 2π * range(-0.25, 0.25, 10001)
const ρ0 = [1 0 0; 0 0 0; 0 0 0]
const ρ1 = [0 0 0; 0 0 0; 0 0 1]
const ρ2 = [0.5 0 0; 0 0 0; 0 0 0.5]

function get_33(det, t)
    sys = get_sys(det)
    return real(propagate(sys, ρ0, t)[3, 3])
end

function get_33_1(det, t)
    sys = get_sys(det)
    return real(propagate(sys, ρ1, t)[3, 3])
end

function get_33_2(det, t)
    sys = get_sys(det)
    return real(propagate(sys, ρ2, t)[3, 3])
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

const prefix = joinpath(@__DIR__, "../imgs/raman-decoherence")

figure()
plot(ts, (t->(real(propagate(sys0, ρ0, t)[1, 1]))).(ts), "C0", label="g1")
plot(ts, (t->(real(propagate(sys0, ρ0, t)[3, 3]))).(ts), "C1", label="g2")
ylim([0, 1])
grid()
legend(fontsize=12)
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_res_t")

figure()
plot(ts, get_33.(2π * 0.1, ts), "C0", label="\$\\delta=2\\pi\\cdot0.1\$")
plot(ts, get_33.(-2π * 0.1, ts), "C2", label="\$\\delta=-2\\pi\\cdot0.1\$")
plot(ts, analytical.(2π * 0.1, ts), "C1", label="\$\\delta=2\\pi\\cdot0.1\$, no decay")
ylim([0, 1])
grid()
legend(fontsize=12)
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_det_t")

figure()
plot(dets ./ (2π), get_33.(dets, 10), "C0", label="t=10")
plot(dets ./ (2π), get_33.(dets, 30), "C1", label="t=30")
plot(dets ./ (2π), get_33.(dets, 70), "C2", label="t=70")
plot(dets ./ (2π), get_33.(dets, 150), "C3", label="t=150")
plot(dets ./ (2π), get_33.(dets, 310), "C4", label="t=310")
plot(dets ./ (2π), get_33.(dets, 630), "C5", label="t=630")
plot(dets ./ (2π), get_33.(dets, 1270), "C6", label="t=1270")
title("\$\\rho_0=|g1\\rangle\\langle g1|\$")
ylim([0, 1])
grid()
legend(fontsize=12, ncol=3)
xlabel("det")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_det1")

figure()
plot(dets ./ (2π), get_33_1.(dets, 10), "C0", label="t=10")
plot(dets ./ (2π), get_33_1.(dets, 30), "C1", label="t=30")
plot(dets ./ (2π), get_33_1.(dets, 70), "C2", label="t=70")
plot(dets ./ (2π), get_33_1.(dets, 150), "C3", label="t=150")
plot(dets ./ (2π), get_33_1.(dets, 310), "C4", label="t=310")
plot(dets ./ (2π), get_33_1.(dets, 630), "C5", label="t=630")
plot(dets ./ (2π), get_33_1.(dets, 1270), "C6", label="t=1270")
title("\$\\rho_0=|g2\\rangle\\langle g2|\$")
ylim([0, 1])
grid()
legend(fontsize=12, ncol=3)
xlabel("det")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_det2")

figure()
plot(dets ./ (2π), get_33_2.(dets, 10), "C0", label="t=10")
plot(dets ./ (2π), get_33_2.(dets, 30), "C1", label="t=30")
plot(dets ./ (2π), get_33_2.(dets, 70), "C2", label="t=70")
plot(dets ./ (2π), get_33_2.(dets, 150), "C3", label="t=150")
plot(dets ./ (2π), get_33_2.(dets, 310), "C4", label="t=310")
plot(dets ./ (2π), get_33_2.(dets, 630), "C5", label="t=630")
plot(dets ./ (2π), get_33_2.(dets, 1270), "C6", label="t=1270")
title("\$\\rho_0=(|g1\\rangle\\langle g1|+|g2\\rangle\\langle g2|)/2\$")
ylim([0, 1])
grid()
legend(fontsize=12, ncol=3)
xlabel("det")
ylabel("p")
NaCsPlot.maybe_save("$(prefix)_det12")

NaCsPlot.maybe_show()

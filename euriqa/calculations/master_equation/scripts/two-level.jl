#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

const sys = System{ComplexF64}(2)
const H = new_H(sys)
H_add_E!(H, 2π * 15, 2)
H_add_Ω!(H, 2π * 10, 1, 2)
add_H!(sys, H)

const ts = range(0, 1, 10001)
const ρ0 = [1 0; 0 0]

function analytical(t)
    Ω = 2π * 10
    Δ = 2π * 15
    Ωg² = Ω^2 + Δ^2
    Ωg = sqrt(Ωg²)
    return Ω^2 / Ωg² * sin(Ωg * t / 2)^2
end

const prefix = joinpath(@__DIR__, "../imgs/two-level")

figure()
plot(ts, (t->(real(propagate(sys, ρ0, t)[2, 2]))).(ts))
plot(ts, analytical.(ts), "--")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

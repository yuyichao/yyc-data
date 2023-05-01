#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

const sys = System{ComplexF64}(2)
const C = new_C(sys)
C_add_decay!(C, 2π * 2.0, 1, 2)
add_C!(sys, C)

const ts = range(0, 1, 10001)
const ρ0 = [0 0; 0 1]

function analytical(t)
    Γ = 2π * 2.0
    return exp(-Γ * t)
end

const prefix = joinpath(@__DIR__, "../imgs/decay")

figure()
plot(ts, (t->(real(propagate(sys, ρ0, t)[1, 1]))).(ts), "C0")
plot(ts, 1 .- analytical.(ts), "C2--")
plot(ts, (t->(real(propagate(sys, ρ0, t)[2, 2]))).(ts), "C1")
plot(ts, analytical.(ts), "C3--")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

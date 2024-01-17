#!/usr/bin/julia

include("util.jl")

using BenchmarkTools
using NaCsPlot
using PyPlot

const params = SystemParams{Float64,0}(Γ=10.0, Ω₁=0.0, Ω₂=0.1, Δ₁=0, Δ₂=0,
                                       ωs=(), ηs=(), nmotions=())
push!(params.decay_branch, DecayParams{Float64,0}(p=1.0, amp=(1, 0), ηs=()))
const builder = Builder(params)

const rates = build_rate_matrices(builder)

const p0 = [0, 0, 1]
const ts = range(0, 3000, 1000)
const result = RateResult(rates, p0)

const p1s = Float64[]
const p2s = Float64[]

for t in ts
    ps = get_probabilities(propagate!(result, t))
    push!(p1s, ps[1])
    push!(p2s, ps[3])
end

plot(ts, p1s, "C0-", label="\$p_1\$")
plot(ts, p2s, "C1-", label="\$p_2\$")
plot(ts, 1 .- exp.(ts .* -0.001), "C2--")
plot(ts, exp.(ts .* -0.001), "C3--")

grid()

show()

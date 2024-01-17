#!/usr/bin/julia

include("util.jl")

using BenchmarkTools
using NaCsPlot
using PyPlot

const params = SystemParams{Float64,0}(Γ=10.0, Ω₁=0.2, Ω₂=0.1, Δ₁=10.0, Δ₂=10.0,
                                       ωs=(), ηs=(), nmotions=())
push!(params.decay_branch, DecayParams{Float64,0}(p=.5, amp=(1, 0), ηs=()))
push!(params.decay_branch, DecayParams{Float64,0}(p=.5, amp=(0, 1), ηs=()))
const builder = Builder(params)

const rates = build_rate_matrices(builder)

const p0 = [0, 0, 1]
const ts = range(0, 6000, 1000)
const result = RateResult(rates, p0)

const p1s = Float64[]
const p2s = Float64[]
const total_rates = Float64[]

for t in ts
    ps = get_probabilities(propagate!(result, t))
    push!(p1s, ps[1])
    push!(p2s, ps[3])
    push!(total_rates, get_total_rate(result))
end

figure()
plot(ts, p1s, "C0-", label="\$p_1\$")
plot(ts, p2s, "C1-", label="\$p_2\$")
legend()
grid()

figure()
plot(ts, total_rates, "C0-")
grid()

show()

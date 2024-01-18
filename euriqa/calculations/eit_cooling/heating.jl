#!/usr/bin/julia

include("util.jl")

using BenchmarkTools
using NaCsPlot
using PyPlot

const nmotion = 400

const params = SystemParams{Float64,1}(Γ=10.0, Ω₁=1.0, Ω₂=0, Δ₁=0.0, Δ₂=0.0,
                                       ωs=(1.0,), ηs=(0.5,), nmotions=(nmotion,))
push!(params.decay_branch, DecayParams{Float64,1}(p=1, amp=(1, 0), ηs=(0.5,)))
const builder = Builder(params)

const rates = @time build_rate_matrices(builder)

const p0 = zeros(3 * nmotion)
p0[1] = 1
const ts = range(0, 10000, 100)
const result = RateResult(rates, p0)

const p1s = Float64[]
const p2s = Float64[]
const pes = Float64[]
const nbars = Float64[]
const nbar_weights = [0:nmotion - 1; 0:nmotion - 1; 0:nmotion - 1]

for t in ts
    ps = get_probabilities(propagate!(result, t))
    push!(p1s, sum(@view ps[1:nmotion]))
    push!(p2s, sum(@view ps[2 * nmotion + 1:3 * nmotion]))
    push!(pes, sum(@view ps[nmotion + 1:2 * nmotion]))
    push!(nbars, dot(ps, nbar_weights))
end

figure()
plot(ts, p1s, "C0-", label="\$p_1\$")
plot(ts, p2s, "C1-", label="\$p_2\$")
plot(ts, pes, "C2-", label="\$p_e\$")
legend()
grid()

figure()
plot(ts, nbars, "C0-")
grid()

show()

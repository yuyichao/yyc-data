#!/usr/bin/julia

using PyPlot

include("utils.jl")

const ψ_2p = zeros(16)
ψ_2p[1] = 0.5
ψ_2p[6] = -0.5
ψ_2p[11] = 0.5
ψ_2p[16] = -0.5
const ψ_1p = zeros(16)
ψ_1p[1] = sqrt(0.5)
ψ_1p[11] = -sqrt(0.5)

function compute_fidelities(εm, εr, εz)
    N = 4
    VN = Val(N)
    O1 = generation_circuit(εr, εz)
    O2 = distillation_circuit(εr, εz)
    ψ1 = O1 * ψ0
    ψ2 = O2 * ψ1
    f1 = fidelity_per_bit(ψ1, ψ_2p)

    ρ0, ρ1 = add_measure_error(project_bit(VN, ψ2, 1)..., εm)
    ρ0, ρ1 = add_measure_error(project_bit(VN, ρ1, 3)..., εm)
    ρ1 = reset_bit(VN, ρ1, 1)
    ρ1 = reset_bit(VN, ρ1, 3)

    p = tr(ρ1)
    ρ1 ./= p
    f2 = fidelity(ρ1, ψ_1p)

    return p, f1, f2
end

const εzs = range(0, 1, 1000)

const ps = Float64[]
const f1s = Float64[]
const f2s = Float64[]

for εz in εzs
    p, f1, f2 = compute_fidelities(0.01, 0.01, εz)
    push!(ps, p)
    push!(f1s, f1)
    push!(f2s, f2)
end

figure()
plot(1 .- f1s, 1 .- f1s)
plot(1 .- f1s, 1 .- f2s)
plot(1 .- f1s, ps)
grid()

figure()
plot(εzs, 1 .- f1s)
plot(εzs, 1 .- f2s)
plot(εzs, ps)
grid()

show()

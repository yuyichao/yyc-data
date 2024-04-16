#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Trap

using ExponentialAction
using LinearAlgebra

function get_x(N, scale)
    M = zeros(ComplexF64, N, N)
    for i in 0:(N - 2)
        M[i + 1, i + 2] = M[i + 2, i + 1] = scale * sqrt(i + 1)
    end
    return M
end

function get_N(N, scale)
    M = zeros(ComplexF64, N, N)
    for i in 1:(N - 1)
        M[i + 1, i + 1] = scale * i
    end
    return M
end

const N = 300

const N_A = 6.02214076e23

const m_Ba138 = 138e-3 / N_A
const λBaP1 = 493e-9
const ωBa138 = 2π * 1e6
const ηBa138P1 = Trap.η(m_Ba138, ωBa138 / (2π), 2π / λBaP1) * sqrt(2)

const m_Yb171 = 171e-3 / N_A
const λYbD3 = 1389e-9
const ωYb171 = 2π * 100e3
const ηYb171D3 = Trap.η(m_Yb171, ωYb171 / (2π), 2π / λYbD3) * sqrt(2)

const iHBa138 = get_N(N, im * ωBa138)
const ikxBa138P1 = get_x(N, im * ηBa138P1)

const iHYb171 = get_N(N, im * ωYb171)
const ikxYb171D3 = get_x(N, im * ηYb171D3)

struct MotionCalculator
    expikx::Matrix{ComplexF64}
    expiHt::Matrix{ComplexF64}
    buff1::Matrix{ComplexF64}
    buff2::Matrix{ComplexF64}
    nbar::Float64
    ω::Float64
    η::Float64
    function MotionCalculator(nbar, ω, η)
        N = ceil(Int, (nbar + 2) * 7) + 1
        return new(exp(get_x(N, im * η)), zeros(ComplexF64, N, N),
                   zeros(ComplexF64, N, N), zeros(ComplexF64, N, N), nbar, ω, η)
    end
end

function get_motion_overlap(calc::MotionCalculator, t)
    expiHt = calc.expiHt
    N = size(expiHt, 1)
    ω = calc.ω
    for i in 1:N
        expiHt[i, i] = cis(ω * t * (i - 1))
    end
    mul!(calc.buff1, expiHt, calc.expikx)
    mul!(calc.buff2, calc.expikx', calc.buff1)
    nbar = calc.nbar
    f = 0.0
    for i in 1:N
        n = i - 1
        p = (nbar / (nbar + 1))^n / (nbar + 1)
        f += p * abs2(calc.buff2[i, i])
    end
    return f
end

function get_fidelity(calc::MotionCalculator, t)
    expiHt = calc.expiHt
    N = size(expiHt, 1)
    ω = calc.ω
    for i in 1:N
        expiHt[i, i] = cis(ω * t * (i - 1))
    end
    mul!(calc.buff1, expiHt, calc.expikx)
    mul!(calc.buff2, calc.expikx', calc.buff1)
    mul!(calc.buff1, calc.expiHt', calc.buff2)
    nbar = calc.nbar
    f = 0.0
    for i in 1:N
        n = i - 1
        p = (nbar / (nbar + 1))^n / (nbar + 1)
        f += p * (1 + real(calc.buff1[i, i])) / 2
    end
    return f
end


struct AnalyticMotionCalculator
    nbar::Float64
    ω::Float64
    η::Float64
    function AnalyticMotionCalculator(nbar, ω, η)
        return new(nbar, ω, η)
    end
end

function get_fidelity(calc::AnalyticMotionCalculator, t)
    η = calc.η
    ω = calc.ω
    nbar = calc.nbar
    η² = η^2
    s, c = sincos(ω * t)
    nd = cis(-η² * s) * exp(-η² * (1 - c) * (2 * nbar + 1))
    return (1 + real(nd)) / 2
end

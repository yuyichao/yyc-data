#!/usr/bin/julia

using QuantumOptics
using LinearAlgebra

struct PulseParam
    Ω::Float64
    δ::Float64
    φ::Float64
    t::Float64
end

struct TwoLevelMotionData{H,S,O} <: Function
    H0::H
    σx::S
    σy::S
    op::O
    pulses::Vector{PulseParam}
    ts::Vector{Float64}
    φ0s::Vector{Float64}

    function TwoLevelMotionData(ω, η, nmax)
        s_basis = SpinBasis(1//2)
        m_basis = FockBasis(nmax)
        basis = s_basis ⊗ m_basis

        h0 = identityoperator(s_basis) ⊗ (number(m_basis) * ω)
        σp = sigmap(s_basis) ⊗ displace(m_basis, im * η)
        σm = sigmam(s_basis) ⊗ displace(m_basis, -im * η)

        σx = σp + σm
        σy = im * σp - im * σm

        op_data = zeros(ComplexF64, (nmax + 1) * 2, (nmax + 1) * 2)
        op = Operator(basis, op_data)

        return new{typeof(h0),typeof(σx),typeof(op)}(h0, σx, σy, op, PulseParam[],
                                                      Float64[0.0], Float64[0.0])
    end
end

function precompute_pulses!(data::TwoLevelMotionData)
    pulses = data.pulses
    npulses = length(pulses)
    ts = data.ts
    φ0s = data.φ0s
    resize!(ts, npulses + 1)
    resize!(φ0s, npulses + 1)
    t = 0.0
    φ0 = 0.0
    for i in 1:npulses
        pulse = pulses[i]
        t += pulse.t
        φ0 += pulse.t * pulse.δ
        ts[i + 1] = t
        φ0s[i + 1] = φ0
    end
end

function (data::TwoLevelMotionData)(t, ψ)
    op = data.op
    op_data = op.data
    op_data .= data.H0.data
    idx = searchsortedlast(data.ts, t)
    npulses = length(data.pulses)
    if idx > npulses
        idx = npulses
    end
    t0 = data.ts[idx]
    φ0 = data.φ0s[idx]
    pulse = data.pulses[idx]
    Ω = pulse.Ω
    if Ω != 0
        θ = pulse.δ * t + pulse.φ + φ0
        s, c = Ω / 2 .* sincos(θ)
        op_data .+= c .* data.σx.data .+ s .* data.σy.data
    end
    return op
end

function run_pulses(data, ψ, npoints=1001; kws...)
    precompute_pulses!(data)
    return timeevolution.schroedinger_dynamic(range(0, data.ts[end], npoints),
                                              ψ, data; kws...)
end

function init_n(data, n)
    s_basis, m_basis = data.op.basis_l.bases
    return spindown(s_basis) ⊗ fockstate(m_basis, n)
end

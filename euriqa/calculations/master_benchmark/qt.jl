#!/usr/bin/julia

module QTTest

using QuantumToolbox
using SparseArrays

rexpects(m, states) = [real(expect(m, s)) for s in states]
function expects(m, states)
    er = Vector{Float64}(undef, length(states))
    ei = Vector{Float64}(undef, length(states))
    for (i, s) in enumerate(states)
        e = expect(m, s)
        er[i] = real(e)
        ei[i] = imag(e)
    end
    return er, ei
end

function _motion_problem(Ω, ω, δ, nbar, nMax)
    η = 2π / 578e-9 * sqrt(1.05e-34 / (2 * 171 * 1.67e-27 * (ω * 1e6)))

    I_s = qeye(2)
    I_m = qeye(nMax)
    n_m = num(nMax)
    ikx_m = im * η * (create(nMax) + destroy(nMax))

    H = ω * I_s ⊗ n_m +
        Ω / 2 * (sigmap() ⊗ exp(-ikx_m) + sigmam() ⊗ exp(ikx_m)) +
        δ * fock_dm(2, 1) ⊗ I_m
    ψ0 = fock_dm(2, 0) ⊗ thermal_dm(nMax, nbar + 1e-100)
    return H, ψ0
end

struct MotionProblem{HT,ΨT}
    H::HT
    ψ0::ΨT
end

function motion_problem(Ω, ω, δ, nbar, nMax)
    return MotionProblem(_motion_problem(Ω, ω, δ, nbar, nMax)...)
end
function motion2_problem(Ω, ω, δ, nbar, nMax)
    H, ψ0 = _motion_problem(Ω, ω, δ, nbar, nMax)
    return MotionProblem(Qobj(sparse(H.data), dims=H.dims), ψ0)
end

function psolve(p::MotionProblem, tlist; kws...)
    return tlist, mesolve(p.H, p.ψ0, tlist; progress_bar=false, kws...).states
end

function pexpect(p::MotionProblem, (tlist, states))
    nMax = size(p.H, 1) ÷ 2
    return (tlist, rexpects(fock_dm(2, 0) ⊗ qeye(nMax), states),
            rexpects(fock_dm(2, 1) ⊗ qeye(nMax), states),
            rexpects(qeye(2) ⊗ num(nMax), states),
            rexpects(qeye(2) ⊗ fock_dm(nMax, 0), states),
            expects(sigmap() ⊗ fock_dm(nMax, 0), states)...)
end

struct SparseProblem{HT,ΨT}
    H::HT
    ψ0::ΨT
    i0::Int
    i1::Int
end

function sparse_problem(H, i0, i1)
    n = size(H, 1)
    return SparseProblem(Qobj(H), fock_dm(n, i0 - 1), i0, i1)
end

function psolve(p::SparseProblem, tlist; kws...)
    return tlist, mesolve(p.H, p.ψ0, tlist; progress_bar=false, kws...).states
end

function pexpect(p::SparseProblem, (tlist, states))
    n = size(p.H, 1)

    return (tlist, rexpects(fock_dm(n, p.i0 - 1), states),
            rexpects(fock_dm(n, p.i1 - 1), states),
            expects(fock(n, p.i0 - 1) * fock(n, p.i1 - 1)', states)...)
end

end

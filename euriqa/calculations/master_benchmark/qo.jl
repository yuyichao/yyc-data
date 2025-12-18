#!/usr/bin/julia

module QOTest

using QuantumOptics
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
    sbasis = SpinBasis(1//2)
    mbasis = FockBasis(nMax - 1)

    I_s = identityoperator(sbasis)
    I_m = identityoperator(mbasis)
    n_m = number(mbasis)
    ikx_m = im * η * (create(mbasis) + destroy(mbasis))

    H = ω * n_m ⊗ I_s +
        Ω / 2 * (exp(-ikx_m) ⊗ sigmap(sbasis) + exp(ikx_m) ⊗ sigmam(sbasis)) +
        δ * I_m ⊗ (spindown(sbasis) ⊗ spindown(sbasis)')
    ψ0 = thermalstate(n_m, nbar + 1e-100) ⊗ (spinup(sbasis) ⊗ spinup(sbasis)')
    return H, ψ0
end

struct MotionProblem{HT,ΨT}
    H::HT
    ψ0::ΨT
end

function motion_problem(Ω, ω, δ, nbar, nMax)
    H, ψ0 = _motion_problem(Ω, ω, δ, nbar, nMax)
    return MotionProblem(Operator(H.basis_l, Matrix(H.data)), ψ0)
end
function motion2_problem(Ω, ω, δ, nbar, nMax)
    H, ψ0 = _motion_problem(Ω, ω, δ, nbar, nMax)
    return MotionProblem(Operator(H.basis_l, sparse(H.data)), ψ0)
end

function psolve(p::MotionProblem, tlist; kws...)
    return timeevolution.master(tlist, p.ψ0, p.H, (); kws...)
end

function pexpect(p::MotionProblem, (tlist, states))
    nMax = size(p.H, 1) ÷ 2
    sbasis = SpinBasis(1//2)
    mbasis = FockBasis(nMax - 1)
    I_s = identityoperator(sbasis)
    I_m = identityoperator(mbasis)
    n_m = number(mbasis)

    return (tlist, rexpects(I_m ⊗ (spinup(sbasis) ⊗ spinup(sbasis)'), states),
            rexpects(I_m ⊗ (spindown(sbasis) ⊗ spindown(sbasis)'), states),
            rexpects(n_m ⊗ I_s, states),
            rexpects((fockstate(mbasis, 0) ⊗ fockstate(mbasis, 0)') ⊗ I_s, states),
            expects((fockstate(mbasis, 0) ⊗ fockstate(mbasis, 0)') ⊗
                (spinup(sbasis) ⊗ spindown(sbasis)'), states)...)
end

struct SparseProblem{HT,ΨT}
    H::HT
    ψ0::ΨT
    i0::Int
    i1::Int
end

function sparse_problem(H, i0, i1)
    n = size(H, 1)
    basis = FockBasis(n - 1)
    return SparseProblem(Operator(basis, H), fockstate(basis, i0 - 1), i0, i1)
end

function psolve(p::SparseProblem, tlist; kws...)
    return timeevolution.master(tlist, p.ψ0, p.H, (); kws...)
end

function pexpect(p::SparseProblem, (tlist, states))
    n = size(p.H, 1)

    basis = FockBasis(n - 1)
    s0 = fockstate(basis, p.i0 - 1)
    s1 = fockstate(basis, p.i1 - 1)

    return (tlist, rexpects(s0 ⊗ s0', states), rexpects(s1 ⊗ s1', states),
            expects(s0 ⊗ s1', states)...)
end

end

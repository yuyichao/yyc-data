#!/usr/bin/julia

module QTTest

using QuantumToolbox

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

function motion_problem(Ω, ω, δ, nbar, nMax)
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

function solve_motion((H, ψ0), tlist; kws...)
    return tlist, mesolve(H, ψ0, tlist; progress_bar=false, kws...).states
end

function expect_motion((H, ψ0), (tlist, states))
    nMax = size(H, 1) ÷ 2

    return (tlist, rexpects(fock_dm(2, 0) ⊗ qeye(nMax), states),
            rexpects(fock_dm(2, 1) ⊗ qeye(nMax), states),
            rexpects(qeye(2) ⊗ num(nMax), states),
            rexpects(qeye(2) ⊗ fock_dm(nMax, 0), states),
            expects(sigmap() ⊗ fock_dm(nMax, 0), states)...)
end

function motion2_problem(Ω, ω, δ, nbar, nMax)
    η = 2π / 578e-9 * sqrt(1.05e-34 / (2 * 171 * 1.67e-27 * (ω * 1e6)))

    I_s = qeye(2)
    I_m = qeye(nMax)
    n_m = num(nMax)
    ikx_m = im * η * (create(nMax) + destroy(nMax))

    H = ω * n_m ⊗ I_s +
        Ω / 2 * (exp(-ikx_m) ⊗ sigmap() + exp(ikx_m) ⊗ sigmam()) +
        δ * I_m ⊗ fock_dm(2, 1)
    ψ0 = thermal_dm(nMax, nbar + 1e-100) ⊗ fock_dm(2, 0)
    return H, ψ0
end

function solve_motion2((H, ψ0), tlist; kws...)
    return tlist, mesolve(H, ψ0, tlist; progress_bar=false, kws...).states
end

function expect_motion2((H, ψ0), (tlist, states))
    nMax = size(H, 1) ÷ 2

    return (tlist, rexpects(qeye(nMax) ⊗ fock_dm(2, 0), states),
            rexpects(qeye(nMax) ⊗ fock_dm(2, 1), states),
            rexpects(num(nMax) ⊗ qeye(2), states),
            rexpects(fock_dm(nMax, 0) ⊗ qeye(2), states),
            expects(fock_dm(nMax, 0) ⊗ sigmap(), states)...)
end

function sparse_problem(H, i0, i1)
    n = size(H, 1)
    return Qobj(H), fock_dm(n, i0 - 1), i0, i1
end

function solve_sparse((H, ψ0, _, _), tlist; kws...)
    return tlist, mesolve(H, ψ0, tlist; progress_bar=false, kws...).states
end

function expect_sparse((H, ψ0, i0, i1), (tlist, states))
    n = size(H, 1)

    return (tlist, rexpects(fock_dm(n, i0 - 1), states),
            rexpects(fock_dm(n, i1 - 1), states),
            expects(fock(n, i0 - 1) * fock(n, i1 - 1)', states)...)
end

end

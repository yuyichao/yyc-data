#!/usr/bin/julia

module QTTest

using QuantumToolbox

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
    meas = [fock_dm(2, 0) ⊗ qeye(nMax), fock_dm(2, 1) ⊗ qeye(nMax),
            qeye(2) ⊗ num(nMax), qeye(2) ⊗ fock_dm(nMax, 0),
            sigmap() ⊗ fock_dm(nMax, 0)]
    result = [expect(mea, states) for mea in meas]
    return (tlist, real.(result[1]), real.(result[2]), real.(result[3]), real.(result[4]),
            real.(result[5]), imag.(result[5]))
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
    meas = [qeye(nMax) ⊗ fock_dm(2, 0), qeye(nMax) ⊗ fock_dm(2, 1),
            num(nMax) ⊗ qeye(2), fock_dm(nMax, 0) ⊗ qeye(2),
            fock_dm(nMax, 0) ⊗ sigmap()]
    result = [expect(mea, states) for mea in meas]
    return (tlist, real.(result[1]), real.(result[2]), real.(result[3]), real.(result[4]),
            real.(result[5]), imag.(result[5]))
end

end

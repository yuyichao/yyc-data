#!/usr/bin/julia

module SimpleTest

using QuantumToolbox
using LinearAlgebra
using SparseArrays

import OrdinaryDiffEq, DiffEqCallbacks

@noinline function dmaster(dρ, H, ρ)
    mul!(dρ, H, ρ, -im, false)
    mul!(dρ, ρ, H, im, true)
    return
end

function solve_master(H, x0, tlist; alg = OrdinaryDiffEq.DP5(), kws...)
    @inline fout(x, t, _) = copy(x)
    tType = float(eltype(tlist))
    out = DiffEqCallbacks.SavedValues(tType, typeof(x0))
    scb = DiffEqCallbacks.SavingCallback(fout, out, saveat=tlist,
                                         tdir = one(eltype(tlist)))

    @inline function dmaster_(dρ, ρ, _, _)
        dmaster(dρ, H, ρ)
        return dρ
    end

    prob = OrdinaryDiffEq.ODEProblem{true}(dmaster_, x0, (convert(tType, tlist[1]),
                                                          convert(tType, tlist[end])))

    OrdinaryDiffEq.solve(prob, alg;
                         reltol = 1.0e-6, abstol = 1.0e-8,
                         save_everystep = false, save_start = false,
                         save_end = false,
                         callback=scb, kws...)
    return out.t, out.saveval
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
    return H.data, ψ0.data
end

function solve_motion((H, ψ0), tlist; kws...)
    return solve_master(H, ψ0, tlist; kws...)
end

function expect_motion((H, ψ0), (tlist, states))
    nMax = size(H, 1) ÷ 2
    meas = [fock_dm(2, 0) ⊗ qeye(nMax), fock_dm(2, 1) ⊗ qeye(nMax),
            qeye(2) ⊗ num(nMax), qeye(2) ⊗ fock_dm(nMax, 0),
            sigmap() ⊗ fock_dm(nMax, 0)]
    result = [[expect(Qobj(mea.data), Qobj(s)) for s in states] for mea in meas]
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
    return H.data, ψ0.data
end

function solve_motion2((H, ψ0), tlist; kws...)
    return solve_master(H, ψ0, tlist; kws...)
end

function expect_motion2((H, ψ0), (tlist, states))
    nMax = size(H, 1) ÷ 2
    meas = [qeye(nMax) ⊗ fock_dm(2, 0), qeye(nMax) ⊗ fock_dm(2, 1),
            num(nMax) ⊗ qeye(2), fock_dm(nMax, 0) ⊗ qeye(2),
            fock_dm(nMax, 0) ⊗ sigmap()]
    result = [[expect(Qobj(mea.data), Qobj(s)) for s in states] for mea in meas]
    return (tlist, real.(result[1]), real.(result[2]), real.(result[3]), real.(result[4]),
            real.(result[5]), imag.(result[5]))
end

end

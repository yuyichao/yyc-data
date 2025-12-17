#!/usr/bin/julia

module QOTest

using QuantumOptics

function motion_problem(Ω, ω, δ, nbar, nMax)
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

function solve_motion((H, ψ0), tlist; kws...)
    return timeevolution.master(tlist, ψ0, H, (); kws...)
end

function motion2_problem(Ω, ω, δ, nbar, nMax)
    η = 2π / 578e-9 * sqrt(1.05e-34 / (2 * 171 * 1.67e-27 * (ω * 1e6)))
    sbasis = SpinBasis(1//2)
    mbasis = FockBasis(nMax - 1)

    I_s = identityoperator(sbasis)
    I_m = identityoperator(mbasis)
    n_m = number(mbasis)
    ikx_m = im * η * (create(mbasis) + destroy(mbasis))

    H = ω * I_s ⊗ n_m +
        Ω / 2 * (sigmap(sbasis) ⊗ exp(-ikx_m) + sigmam(sbasis) ⊗ exp(ikx_m)) +
        δ * (spindown(sbasis) ⊗ spindown(sbasis)') ⊗ I_m
    ψ0 = (spinup(sbasis) ⊗ spinup(sbasis)') ⊗ thermalstate(n_m, nbar + 1e-100)
    return H, ψ0
end

function solve_motion2((H, ψ0), tlist; kws...)
    return timeevolution.master(tlist, ψ0, H, (); kws...)
end

end

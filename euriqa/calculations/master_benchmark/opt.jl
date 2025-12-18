#!/usr/bin/julia

module OPTTest

using QuantumToolbox
using LinearAlgebra
using SparseArrays

import OrdinaryDiffEq, DiffEqCallbacks

@noinline function dmaster(dρ, H, ρ)
    mul!(dρ, H, ρ, -im, false)
    mul!(dρ, ρ, H, im, true)
    return
end

@inline function _apply_hc!(dρ)
    n = size(dρ, 1)
    # compute -i * dρ^dagger + i * dρ
    @inbounds for i in 1:n
        dρ[i + (i - 1) * n] = -2 * imag(dρ[i + (i - 1) * n])
        for j in (i + 1):n
            idx1 = j + (i - 1) * n
            idx2 = i + (j - 1) * n
            v1 = dρ[idx1]
            v2 = dρ[idx2]

            re_o = -imag(v1) - imag(v2)
            im_o = real(v2) - real(v1)

            dρ[idx1] = complex(re_o, -im_o)
            dρ[idx2] = complex(re_o, im_o)
        end
    end
end

@noinline function dmaster(dρ::DenseMatrix, H::DenseMatrix, ρ::DenseMatrix)
    mul!(dρ, ρ, H, true, false)
    _apply_hc!(dρ)
end
@noinline function dmaster(dρ::DenseMatrix, H::SparseArrays.AbstractSparseMatrixCSC, ρ::DenseMatrix)
    if nnz(H) < size(dρ, 1)
        # If H is really sparse, the matrix multiplication could be cheaper than
        # computing the hermitian conjugate.
        mul!(dρ, H, ρ, -im, false)
        mul!(dρ, ρ, H, im, true)
    else
        mul!(dρ, ρ, H, true, false)
        _apply_hc!(dρ)
    end
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

function rexpects(m, states)
    m = Qobj(m.data)
    [real(expect(m, Qobj(s))) for s in states]
end
function expects(m, states)
    m = Qobj(m.data)
    er = Vector{Float64}(undef, length(states))
    ei = Vector{Float64}(undef, length(states))
    for (i, s) in enumerate(states)
        e = expect(m, Qobj(s))
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
    return H.data, ψ0.data
end

function solve_motion((H, ψ0), tlist; kws...)
    return solve_master(H, ψ0, tlist; kws...)
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
    return H.data, ψ0.data
end

function solve_motion2((H, ψ0), tlist; kws...)
    return solve_master(H, ψ0, tlist; kws...)
end

function expect_motion2((H, ψ0), (tlist, states))
    nMax = size(H, 1) ÷ 2

    return (tlist, rexpects(qeye(nMax) ⊗ fock_dm(2, 0), states),
            rexpects(qeye(nMax) ⊗ fock_dm(2, 1), states),
            rexpects(num(nMax) ⊗ qeye(2), states),
            rexpects(fock_dm(nMax, 0) ⊗ qeye(2), states),
            expects(fock_dm(nMax, 0) ⊗ sigmap(), states)...)
end

end

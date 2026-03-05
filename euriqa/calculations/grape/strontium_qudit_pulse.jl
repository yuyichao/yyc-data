#!/usr/bin/julia

include("grape.jl")

using JuMP
using NLopt

using AMO
using AMO: Atomic

using LinearAlgebra
using StaticArrays

const OPType = Matrix{ComplexF64}

const μB = 1.39962449171
const μN = 0.7622593285e-3

const g_I = -1.093 * μN / (9/2) / μB
const g_3P1 = AMO.g_sum(1, 1, AMO.g_l, 1, AMO.g_s)

ground_matrix(B) = Atomic.hyperfine_matrix(I=9/2, J=0, Bm=B, g_I=2π * g_I * μB, g_J=0,
                                           Ahf=0, Bhf=0)
excited_matrix(B) = Atomic.hyperfine_matrix(I=9/2, J=1, Bm=B, g_I=2π * g_I * μB,
                                            g_J=2π * g_3P1 * μB,
                                            Ahf=2π * -260.083, Bhf=2π * -35.355)
couple_matrix(Ωs) = Atomic.dipole_couple_matrix(0, 1, Ωs; S=9/2)

mutable struct Drive <: AbstractOP{OPType,2}
    base::Matrix{ComplexF64}
    amp_term::Matrix{ComplexF64}
    det_term::Matrix{ComplexF64}

    buff::Matrix{ComplexF64}
    buff2::Matrix{ComplexF64}

    grad_amp::Matrix{ComplexF64}
    grad_det::Matrix{ComplexF64}

    amp::Float64
    det::Float64

    function Drive(pol, B)
        G = ground_matrix(B)
        E = excited_matrix(B)
        C = couple_matrix(pol)

        sz_g = size(G, 1)
        sz_e = size(E, 1)
        sz = sz_g + sz_e
        @assert size(G) == (sz_g, sz_g)
        @assert size(E) == (sz_e, sz_e)
        @assert size(C) == (sz_g, sz_e)

        base = zeros(ComplexF64, sz, sz)
        amp = zeros(ComplexF64, sz, sz)
        det = zeros(ComplexF64, sz, sz)

        base[1:sz_g, 1:sz_g] .= G
        base[sz_g + 1:sz_g + sz_e, sz_g + 1:sz_g + sz_e] .= E

        amp[1:sz_g, sz_g + 1:sz_g + sz_e] .= C
        amp[sz_g + 1:sz_g + sz_e, 1:sz_g] .= C'

        for i in sz_g + 1:sz_g + sz_e
            det[i, i] = 1
        end
        return Kernel(base, amp, det, zeros(ComplexF64, sz, sz),
                      zeros(ComplexF64, sz * 2, sz * 2), 0, 0)
    end
end
function set_params(dri::Drive, params)
    dri.amp = params[1]
    dri.det = params[2]
    return
end

function compute(dri::Drive, grad)
    amp = dri.amp
    det = dri.det

    sz = size(dri.base, 1)
    buff = dri.buff
    buff .= dri.base .+ amp .* dri.amp_term .+ det .* dri.det_term

    buff2 = dri.buff2

    buff2[1:sz, 1:sz] .= buff
    buff2[sz + 1:sz * 2, sz + 1:sz * 2] .= buff
    buff2[1:sz, sz + 1:sz * 2] .= dri.amp_term
    buff2 = LinearAlgebra.exp!(buff2)
    dri.grad_amp .= @view(buff2[sz + 1:sz * 2, sz + 1:sz * 2])

    buff2[1:sz, 1:sz] .= buff
    buff2[sz + 1:sz * 2, sz + 1:sz * 2] .= buff
    buff2[1:sz, sz + 1:sz * 2] .= dri.det_term
    buff2 = LinearAlgebra.exp!(buff2)
    dri.grad_det .= @view(buff2[sz + 1:sz * 2, sz + 1:sz * 2])

    buff .= @view(buff2[1:sz, 1:sz])

    grad[1] = dri.grad_amp
    grad[2] = dri.grad_det

    return buff
end

@inline function convert_res(op)
    res = 0.0
    for i in 1:10, j in 1:10
        if i in (1, 2) && j in (1, 2)
            exp_val2 = 0.5
        elseif i == j
            exp_val2 = 1.0
        else
            exp_val2 = 0.0
        end
        res += abs2(abs2(op[i, j]) - exp_val2)
    end
    return res
end
@inline function convert_grad(op, op_grad)
    res = 0.0
    for i in 1:10, j in 1:10
        if i in (1, 2) && j in (1, 2)
            exp_val2 = 0.5
        elseif i == j
            exp_val2 = 1.0
        else
            exp_val2 = 0.0
        end
        v = op[i, j]
        g = op_grad[i, j]
        diff = abs2(v) - exp_val2
        diff_g = 2 * (real(v) * real(g) + imag(v) * imag(g))
        res += 2 * diff * diff_g
    end
    return res
end

@inline function convert_res_grads(op, op_grads, grads)
    grads .= convert_grad.(Ref(op), op_grads)
    return convert_res(op)
end

mutable struct QuditSeq{N,S,P,OB}
    const s::S # Sequence
    const params::P # Input parameters
    res::Float64 # Output result
    const grads::P # Output gradients

    const op_buff::OB

    function QuditSeq{N}(pol, B) where N
        ops = ntuple(_->Drive(pol, B), Val(N + 1))
        s = Sequence{OPType}(ops)
        params = MVector{2N,Float64}(undef)
        grads = MVector{2N,Float64}(undef)
        op_buff = Vector{OPType}(undef, 2N)
        return new{N,typeof(s),typeof(params),typeof(op_buff)}(
            s, params, NaN, grads, op_buff)
    end
end

function update_params!(ps::QuditSeq{N}, params) where N
    @assert length(params) == 2N
    params_ary = ps.params
    has_diff = false
    @inbounds @simd ivdep for i in 1:2N
        p0 = params_ary[i]
        p1 = params[i]
        has_diff |= p0 != p1
        params_ary[i] = p1
    end
    if !isnan(ps.res) && !has_diff
        return
    end
    ps.res = NaN
    set_params(ps.s, ps.params)
    op = compute(ps.s, ps.op_buff)
    ps.res = convert_res_grads(op, ps.op_buff, ps.grads)
    return
end

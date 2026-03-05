#!/usr/bin/julia

using AMO
using AMO: Atomic, TimeSequence as TS

using LinearAlgebra
using StaticArrays

using JuMP
using NLopt

const OPType = Matrix{ComplexF64}

const μB = 1.39962449171
const μN = 0.7622593285e-3

const g_I = -1.093 * μN / (9/2) / μB
const g_3P1 = AMO.g_sum(1, 1, AMO.g_l, 1, AMO.g_s)

ground_matrix(B) = Atomic.hyperfine_matrix(Float64,
                                           I=9/2, J=0, Bm=B, g_I=2π * g_I * μB, g_J=0,
                                           Ahf=0, Bhf=0)
excited_matrix(B) = Atomic.hyperfine_matrix(Float64,
                                            I=9/2, J=1, Bm=B, g_I=2π * g_I * μB,
                                            g_J=2π * g_3P1 * μB,
                                            Ahf=2π * -260.083, Bhf=2π * -35.355)
couple_matrix(Ωs) = Atomic.dipole_couple_matrix(0, 1, Ωs; S=9/2)

mutable struct Kernel
    const base::Matrix{ComplexF64}
    const amp_term::Matrix{ComplexF64}

    const iH_buff::Matrix{ComplexF64}
    const iHt_buff::Matrix{ComplexF64}
    const buff2::Matrix{ComplexF64}

    const sz_g::Int
    const sz_e::Int

    function Kernel(pol, B)
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

        base[1:sz_g, 1:sz_g] .= complex.(0, G)
        base[sz_g + 1:sz_g + sz_e, sz_g + 1:sz_g + sz_e] .= complex.(0, E)

        amp[1:sz_g, sz_g + 1:sz_g + sz_e] .= im * C
        amp[sz_g + 1:sz_g + sz_e, 1:sz_g] .= im * C'
        return new(base, amp, zeros(ComplexF64, sz, sz),
                   zeros(ComplexF64, sz, sz), zeros(ComplexF64, sz * 2, sz * 2),
                   sz_g, sz_e)
    end
end

function compute_values!(kern::Kernel, amp, det, t, res, grad_amp, grad_det, grad_t)
    sz_g = kern.sz_g
    sz_e = kern.sz_e
    sz = sz_g + sz_e

    base = kern.base
    amp_term = kern.amp_term

    iH = kern.iH_buff
    @inbounds for j in 1:sz_g
        for i in 1:sz_g
            iH[i, j] = base[i, j]
        end
        for i in sz_g + 1:sz
            iH[i, j] = amp_term[i, j] * amp
        end
    end
    @inbounds for j in sz_g + 1:sz
        for i in 1:sz_g
            iH[i, j] = amp_term[i, j] * amp
        end
        for i in sz_g + 1:sz
            iH[i, j] = base[i, j]
        end
        iH[j, j] += complex(0, det)
    end

    iHt = kern.iHt_buff
    iHt .= iH .* t

    buff2 = kern.buff2

    buff2[1:sz, 1:sz] .= iHt
    buff2[sz + 1:sz * 2, sz + 1:sz * 2] .= iHt
    buff2[1:sz, sz + 1:sz * 2] .= amp_term .* t
    buff2 = LinearAlgebra.exp!(buff2)
    grad_amp .= @view(buff2[1:sz, sz + 1:sz * 2])

    buff2[1:sz, 1:sz] .= iHt
    buff2[sz + 1:sz * 2, sz + 1:sz * 2] .= iHt
    buff2[1:sz, sz + 1:sz * 2] .= 0
    @inbounds for j in sz_g + 1:sz
        buff2[j, sz + j] = complex(0, t)
    end
    buff2 = LinearAlgebra.exp!(buff2)
    grad_det .= @view(buff2[1:sz, sz + 1:sz * 2])

    res .= @view(buff2[1:sz, 1:sz])
    mul!(grad_t, iH, res)
    return
end

mutable struct Drive <: TS.AbstractStep{OPType,3}
    const kern::Kernel
    const res::Matrix{ComplexF64}

    const grad_amp::Matrix{ComplexF64}
    const grad_det::Matrix{ComplexF64}
    const grad_t::Matrix{ComplexF64}

    amp::Float64
    det::Float64
    t::Float64

    function Drive(kern::Kernel)
        sz = kern.sz_g + kern.sz_e
        return new(kern, zeros(ComplexF64, sz, sz), zeros(ComplexF64, sz, sz),
                   zeros(ComplexF64, sz, sz), zeros(ComplexF64, sz, sz), 0, 0, 0)
    end
end

# TS.support_inplace_compute(::Type{Drive}) = true
function TS.set_params!(dri::Drive, params)
    dri.amp = params[1]
    dri.det = params[2]
    dri.t = params[3]
    return
end
function TS.compute(dri::Drive, grads)
    if !isempty(grads)
        grads[1] = dri.grad_amp
        grads[2] = dri.grad_det
        grads[3] = dri.grad_t
    end
    compute_values!(dri.kern, dri.amp, dri.det, dri.t,
                    dri.res, dri.grad_amp, dri.grad_det, dri.grad_t)
    return dri.res
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
        kern = Kernel(pol, B)
        ops = ntuple(_->Drive(kern), Val(N))
        s = TS.Sequence{OPType}(ops)
        params = MVector{3N,Float64}(undef)
        grads = MVector{3N,Float64}(undef)
        op_buff = Vector{OPType}(undef, 3N)
        return new{N,typeof(s),typeof(params),typeof(op_buff)}(
            s, params, NaN, grads, op_buff)
    end
end

function update_params!(ps::QuditSeq{N}, params) where N
    @assert length(params) == 3N
    params_ary = ps.params
    has_diff = false
    @inbounds @simd ivdep for i in 1:3N
        p0 = params_ary[i]
        p1 = params[i]
        has_diff |= p0 != p1
        params_ary[i] = p1
    end
    if !isnan(ps.res) && !has_diff
        return
    end
    ps.res = NaN
    TS.set_params!(ps.s, ps.params)
    op = TS.compute(ps.s, ps.op_buff)
    ps.res = convert_res_grads(op, ps.op_buff, ps.grads)
    return
end

mutable struct SeqModel{M,S}
    const m::M
    const sys::S
    const amp_args::Vector{VariableRef}
    const det_args::Vector{VariableRef}
    const t_args::Vector{VariableRef}
    res::Any
    obj::Any
    function SeqModel(ps::QuditSeq{N}; m=nothing) where N
        if m === nothing
            m = Model(NLopt.Optimizer)
            set_attribute(m, "algorithm", :LD_SLSQP)
        end
        function res_func(params...)
            update_params!(ps, params)
            return ps.res
        end
        function grad_func(g, params...)
            update_params!(ps, params)
            g .= ps.grads
            return
        end
        @variable(m, amp[i=1:N])
        @variable(m, det[i=1:N])
        @variable(m, t[i=1:N])
        register(m, :fidelity, 3N, res_func, grad_func, autodiff=false)
        return new{typeof(m),typeof(ps)}(m, ps, amp, det, t, nothing, nothing)
    end
end

function finalize!(model::SeqModel)
    m = model.m
    N = length(model.amp_args)
    args = vec([v[i] for v in (model.amp_args, model.det_args, model.t_args), i = 1:N])
    res = @NLexpression(m, fidelity(args...))
    @NLobjective(m, Min, res)
    model.res = res
    model.obj = res
    return
end

function optimize_pulse!(model::SeqModel, init_amps, init_dets, init_ts)
    for (var, val) in zip(model.amp_args, init_amps)
        set_start_value(var, val)
    end
    for (var, val) in zip(model.det_args, init_dets)
        set_start_value(var, val)
    end
    for (var, val) in zip(model.t_args, init_ts)
        set_start_value(var, val)
    end
    JuMP.optimize!(model.m)
    return (value(model.obj), value(model.res),
            value.(model.amp_args), value.(model.det_args), value.(model.t_args))
end

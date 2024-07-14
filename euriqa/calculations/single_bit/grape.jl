#!/usr/bin/julia

using StaticArrays

abstract type AbstractOP end

function get_nparams end
function get_op_grads end

struct XYRotation <: AbstractOP
end

get_nparams(op::XYRotation) = 2
function get_op_grads(op::XYRotation, (θ, φ)::NTuple{2})
    sinθ, cosθ = sincos(θ)
    sinφ_2, cosφ_2 = sincos(φ / 2)
    U = SA[cosφ_2 sinφ_2 * (im * cosθ + sinθ)
           sinφ_2 * (im * cosθ - sinθ) cosφ_2]
    gθ = SA[0 sinφ_2 * (-im * sinθ + cosθ)
             sinφ_2 * (-im * sinθ - cosθ) 0]
    gφ = SA[-sinφ_2 / 2 cosφ_2 / 2 * (im * cosθ + sinθ)
             cosφ_2 / 2 * (im * cosθ - sinθ) -sinφ_2 / 2]
    return U, gθ, gφ
end

struct ZRotation <: AbstractOP
end

get_nparams(op::ZRotation) = 1
function get_op_grads(op::ZRotation, (θ,)::NTuple{1})
    sinθ_2, cosθ_2 = sincos(θ / 2)
    U = SA[cosθ_2 + im * sinθ_2 0
           0 cosθ_2 - im * sinθ_2]
    gθ = SA[(-sinθ_2 + im * cosθ_2) / 2 0
             0 (-sinθ_2 - im * cosθ_2) / 2]
    return U, gθ
end

struct GlobalRotation <: AbstractOP
end

get_nparams(op::GlobalRotation) = 1
function get_op_grads(op::GlobalRotation, (θ,)::NTuple{1})
    sinθ, cosθ = sincos(θ)
    U = SA[cosθ + im * sinθ 0
           0 cosθ + im * sinθ]
    gθ = SA[-sinθ + im * cosθ 0
             0 -sinθ + im * cosθ]
    return U, gθ
end

struct CompositeOP{T<:Tuple{Vararg{AbstractOP}}} <: AbstractOP
    ops::T
end

get_nparams(op::CompositeOP) = sum(get_nparams.(op.ops))

@inline function _get_op_grads(ops::NTuple{N,AbstractOP}, params::Tuple) where N
    if N == 1
        return get_op_grads(ops[1], params)
    end
    N1 = N ÷ 2
    ops1 = ops[1:N1]
    ops2 = ops[N1 + 1:N]
    nparams1 = sum(get_nparams.(ops1))
    params1 = params[1:nparams1]
    params2 = params[nparams1 + 1:end]

    op1, grads1... = _get_op_grads(ops1, params1)
    op2, grads2... = _get_op_grads(ops2, params2)

    nparams1 = length(grads1)
    nparams2 = length(grads2)

    let op1 = op1, op2 = op2
        return (op2 * op1, ntuple(@inline(i->op2 * grads1[i]), Val(nparams1))...,
                ntuple(@inline(i->grads2[i] * op1), Val(nparams2))...)
    end
end

function get_op_grads(op::CompositeOP, params::Tuple)
    return _get_op_grads(op.ops, params)
end

mutable struct SingleBitOptimizer{NXY,OP,U,S,NParam}
    const op::OP
    const tgt::U
    const scales::S
    params::NTuple{NParam,Float64}
    diff::Float64
    diff_grad::NTuple{NParam,Float64}
    @inline function SingleBitOptimizer(tgt::U, s::S, nxy) where {U, S}
        ops = (ntuple(i->XYRotation(), nxy)..., ZRotation(), GlobalRotation())
        op = CompositeOP(ops)
        nparams = get_nparams(op)
        opt = new{nxy,typeof(op),U,S,nparams}(op, tgt, s)
        _update_params!(opt, ntuple(i->0.0, Val(nparams)))
        return opt
    end
end

function _update_params!(opt::SingleBitOptimizer{NXY}, params) where NXY
    nparams = length(params)
    diff_total = 0.0
    diff_total_grad = ntuple(i->0.0, Val(nparams))
    for (s, w) in opt.scales
        params_s = ntuple(i->(i % 2 == 0 && 2 * i <= NXY ? params[i] * s : params[i]),
                          Val(nparams))
        op, grads... = get_op_grads(opt.op, params_s)
        diff_u = op .- opt.tgt
        diff = sum(abs2.(diff_u))
        diff_grad = ntuple(i->(2 * sum(real.(diff_u .* conj.(grads[i])))), Val(nparams))

        diff_total += diff * w
        diff_total_grad = diff_total_grad .+ diff_grad .* w
    end
    opt.diff = diff_total
    opt.diff_grad = diff_total_grad
    opt.params = params
    return
end

function update_params!(opt::SingleBitOptimizer, params)
    if params != opt.params
        _update_params!(opt, params)
    end
    return
end

function objective_function(opt::SingleBitOptimizer, params)
    update_params!(opt, params)
    return opt.diff
end
function gradient_function(g, opt::SingleBitOptimizer, params)
    update_params!(opt, params)
    g .= opt.diff_grad
    return
end

using JuMP

function prep_model(opt::SingleBitOptimizer{NXY,OP,U,S,NParam}, model::Model)
    @variable(model, params[i=1:NParam])
    for i in 1:NXY
        set_lower_bound(params[i * 2 - 1], -2π)
        set_upper_bound(params[i * 2 - 1], 2π)
        set_lower_bound(params[i * 2], 0)
        set_upper_bound(params[i * 2], 4π)
    end
    set_lower_bound(params[NParam - 1], -2π)
    set_upper_bound(params[NParam - 1], 2π)
    set_lower_bound(params[NParam], -2π)
    set_upper_bound(params[NParam], 2π)
    register(model, :f, NParam, (params...)->objective_function(opt, params),
             (g, params...)->gradient_function(g, opt, params), autodiff=false)
    @NLobjective(model, Min, f(params...))
    return params
end

function optimize!(opt::SingleBitOptimizer{NXY,OP,U,S,NParam},
                   model::Model, init_params) where {NXY,OP,U,S,NParam}
    params = prep_model(opt, model)
    for i in 1:length(init_params)
        set_start_value(params[i], init_params[i])
    end
    JuMP.optimize!(model)
    paramsv = value.(params)
    objv = objective_function(opt, (paramsv...,))
    return paramsv, objv
end


# 0.7674911777136161, 2.080632852404316, -2.296502015747654, 3.8161185539170317, 0.7541557777311648, 1.7688717934276181, -0.24249358600417162, 2.089992998295301, 0.3931511747957044, 2.1961638351245236, 0.46868099008124475, 2.1025367442608673, -0.5481025256441294, 2.2110297327752075, -0.47880990906321685, 2.227593109649211, 0.8211884262082709, 2.12078523130119, 0.24427209124677646, 2.2028227284270367, -0.8734985619664251, 2.285648961478061, -0.17984388156499567, 2.1399302268456935, 0.5457587390860794, 1.5707937921104242

# -5.770399283358007e-9, 0.10461513164859655, -1.3219957126760782e-8, 0.18101817956823482, 4.301715554723471e-8, 4.817396643675329, -7.804101989147135e-8, 5.6340907050867175, 2.721913812760922e-8, 5.634090682225868, 4.3349802724621556e-8, 5.6340907789016335, 1.1875445052316024e-8, 1.570796366127592

# 0.6572752107353826, 0.4078672060419158, -0.7243747592914618, 8.583243950289846, -3.6984992698108847, 8.477635043411707, -2.6685244304514706, 9.691404679775392, 0.6512655921051046, -1.5708197652011708

# 1.3580893689681928, 0.39295281572760654, -0.37249591866182835, 0.47534812910219615, -5.698312330252018, 10.75425516029681, -2.4267022208604683, 8.127329033646664, 1.2204460747825752, -1.570796145487943

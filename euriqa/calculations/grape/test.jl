#!/usr/bin/julia

include("grape.jl")

mutable struct OP1{VT} <: AbstractOP{VT,1}
    v1::Float64
    OP1{VT}() where VT = new()
end
function set_params(op::OP1, params::NTuple{1})
    op.v1 = params[1]
    return
end
function compute(op::OP1, grad)
    grad[1] = 1
    return op.v1
end
mutable struct OP2{VT} <: AbstractOP{VT,2}
    v1::Float64
    v2::Float64
    OP2{VT}() where VT = new()
end
function set_params(op::OP2, params::NTuple{2})
    op.v1 = params[1]
    op.v2 = params[2]
    return
end
function compute(op::OP2, grad)
    grad[1] = op.v2
    grad[2] = op.v1
    return op.v1 * op.v2
end
mutable struct OP3{VT} <: AbstractOP{VT,3}
    v1::Float64
    v2::Float64
    v3::Float64
    OP3{VT}() where VT = new()
end
function set_params(op::OP3, params::NTuple{3})
    op.v1 = params[1]
    op.v2 = params[2]
    op.v3 = params[3]
    return
end
function compute(op::OP3, grad)
    grad[1] = op.v2 * op.v3
    grad[2] = op.v1 * op.v3
    grad[3] = op.v1 * op.v2
    return op.v1 * op.v2 * op.v3
end

ops = (OP1{Float64}(), OP2{Float64}(), OP3{Float64}(), OP1{Float64}(),
       OP1{Float64}(), OP2{Float64}(), OP3{Float64}(), OP2{Float64}(),
       OP2{Float64}(), OP1{Float64}(), OP3{Float64}(), OP3{Float64}())
s = Sequence{Float64}(ops)
NP = nparams(typeof(s))
params = [1.0 for i in 1:NP]
set_params(s, (params...,))
grads = zeros(NP)
@show compute(s, grads)
@show grads

params[1] = 2.0

set_params(s, (params...,))
@show compute(s, grads)
@show grads

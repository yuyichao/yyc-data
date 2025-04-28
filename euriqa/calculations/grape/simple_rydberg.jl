#/usr/bin/julia

include("grape.jl")

using StaticArrays

const OPType = SMatrix{4,4,ComplexF64}

mutable struct Drive <: AbstractOP{OPType,2}
    θ::Float64
    ϕ::Float64
    Drive() = new()
end
function set_params(dri::Drive, params)
    dri.θ = params[1]
    dri.ϕ = params[2]
    return
end
function compute(dri::Drive, grad)
    st, ct = sincos(dri.θ)
    st2, ct2 = sincos(√2 * dri.θ)

    sp, cp = sincos(dri.ϕ)

    diag1 = ct
    diag1_θ = -st
    diag1_ϕ = 0
    off1 = st * complex(sp, -cp)
    off1_θ = ct * complex(sp, -cp)
    off1_ϕ = st * complex(cp, sp)

    diag2 = ct2
    diag2_θ = -√2 * st2
    diag2_ϕ = 0
    off2 = st2 * complex(sp, -cp)
    off2_θ = √2 * ct2 * complex(sp, -cp)
    off2_ϕ = st2 * complex(cp, sp)

    grad[1] = @SMatrix[diag1_θ       off1_θ   0             0
                       conj(off1_θ)  diag1_θ  0             0
                       0             0        diag2_θ       off2_θ
                       0             0        conj(off2_θ)  diag2_θ]
    grad[2] = @SMatrix[diag1_ϕ       off1_ϕ   0             0
                       conj(off1_ϕ)  diag1_ϕ  0             0
                       0             0        diag2_ϕ       off2_ϕ
                       0             0        conj(off2_ϕ)  diag2_ϕ]
    return @SMatrix[diag1       off1   0           0
                    conj(off1)  diag1  0           0
                    0           0      diag2       off2
                    0           0      conj(off2)  diag2]
end

mutable struct PhaseZ <: AbstractOP{OPType,1}
    ϕ::Float64
    PhaseZ() = new()
end
function set_params(p::PhaseZ, params)
    p.ϕ = params[1]
    return
end
function compute(p::PhaseZ, grad)
    s, c = sincos(p.ϕ)
    s2 = 2 * s * c
    c2 = 2 * c * c - 1

    diag1 = complex(c, s)
    diag1_ϕ = complex(-s, c)

    diag2 = complex(c2, s2)
    diag2_ϕ = 2 * complex(-s2, c2)

    grad[1] = @SMatrix[diag1_ϕ  0  0        0
                       0        0  0        0
                       0        0  diag2_ϕ  0
                       0        0  0        0]
    return @SMatrix[diag1  0  0      0
                    0      1  0      0
                    0      0  diag2  0
                    0      0  0      1]
end

function convert_res(op)
    return abs2(op[1, 1] - 1) * 2 + abs2(op[3, 3] + 1)
end
function convert_grad(op, op_grad)
    op1 = op[1, 1]
    op3 = op[3, 3]
    op1_grad = op_grad[1, 1]
    op3_grad = op_grad[3, 3]

    return (4 * ((real(op1) - 1) * real(op1_grad) + imag(op1) * imag(op1_grad)) +
        2 * ((real(op3) + 1) * real(op3_grad) + imag(op3) * imag(op3_grad)))
end

function convert_res_grads(op, op_grads, grads)
    grads .= convert_grad.(Ref(op), op_grads)
    return convert_res(op)
end

mutable struct PMPulseSeq{N,S,P,OB,RB}
    const s::S # Sequence
    const params::P # Input parameters
    res::Float64 # Output result
    const grads::P # Output gradients

    const op_buff::OB
    const res_buff::RB

    function PMPulseSeq{N}() where N
        ops = ntuple(Val(N + 1)) do i
            return i == 1 ? PhaseZ() : Drive()
        end
        s = Sequence{OPType}(ops)
        params = MVector{N + 2,Float64}(undef)
        grads = MVector{N + 2,Float64}(undef)
        op_buff = Vector{OPType}(undef,2N + 1)
        res_buff = MVector{2N + 1,Float64}(undef)
        return new{N,typeof(s),typeof(params),typeof(op_buff),typeof(res_buff)}(
            s, params, NaN, grads, op_buff, res_buff)
    end
end

function update_params!(ps::PMPulseSeq{N}, params) where N
    if !isnan(ps.res) && all(ps.params .== params)
        return
    end
    res_buff = ps.res_buff
    res_buff[1] = params[1]
    angle = params[2] / N
    for i in 1:N
        res_buff[2 * i] = angle
        res_buff[2 * i + 1] = params[i + 2]
    end

    set_params(ps.s, res_buff)
    op = compute(ps.s, ps.op_buff)
    res = convert_res_grads(op, ps.op_buff, res_buff)

    grads = ps.grads
    grads[1] = res_buff[1]
    angle_grad = 0.0
    for i in 1:N
        angle_grad += res_buff[2 * i]
        grads[i + 2] = res_buff[2 * i + 1]
    end
    grads[2] = angle_grad / N
    ps.res = res
    ps.params .= params
    return
end

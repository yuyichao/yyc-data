#/usr/bin/julia

const OPType = SMatrix{4,4,ComplexF64}

include("grape.jl")

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

    diag1 = c
    diag1_θ = -s
    diag1_ϕ = 0
    off1 = s * complex(sp, -cp)
    off1_θ = c * complex(sp, -cp)
    off1_ϕ = s * complex(cp, sp)

    diag2 = c2
    diag2_θ = -√2 * s2
    diag2_ϕ = 0
    off2 = s2 * complex(sp, -cp)
    off2_θ = √2 * c2 * complex(sp, -cp)
    off2_ϕ = s2 * complex(cp, sp)

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

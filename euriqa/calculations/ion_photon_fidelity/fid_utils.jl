#!/usr/bin/julia

using ForwardDiff
using QuantumOptics
using LinearAlgebra
using Markdown

struct System{N,S1,S2,P1,P2}
    basis1::NLevelBasis{Int}
    ψ1s::Vector{S1}
    ψ2s::Vector{S2}
    proj1s::Vector{P1}
    proj1comps::Vector{P1}
    proj2s::Matrix{P2}
    proj2sub2::Matrix{P2}
    proj2oneof::Vector{P2}
    function System{N}() where N
        basis1 = NLevelBasis(N)
        ψ1s = [nlevelstate(basis1, i) for i in 1:N]
        S1 = eltype(ψ1s)
        ψ2s = [ψ1 ⊗ ψ1 for ψ1 in ψ1s]
        S2 = eltype(ψ2s)
        proj1s = [projector(ψ1) for ψ1 in ψ1s]
        P1 = eltype(proj1s)
        proj1full = sum(proj1s)
        proj1comps = [proj1full - proj for proj in proj1s]
        proj2s = [proj1 ⊗ proj2 for proj1 in proj1s, proj2 in proj1s]
        P2 = eltype(proj2s)
        proj2sub2 = [proj2s[i, i] + proj2s[i, j] + proj2s[j, i] + proj2s[j, j]
                     for i in 1:N, j in 1:N]
        proj2oneof = [proj1s[i] ⊗ proj1comps[i] + proj1comps[i] ⊗ proj1s[i] for i in 1:N]
        return new{N,S1,S2,P1,P2}(basis1, ψ1s, ψ2s, proj1s, proj1comps,
                                  proj2s, proj2sub2, proj2oneof)
    end
end

@inline function Rele(ia, ib, i, j, θ)
    if (i == ia && j == ia) || (i == ib && j == ib)
        return complex(cos(θ / 2))
    elseif (i == ia && j == ib) || (i == ib && j == ia)
        return complex(0, -sin(θ / 2))
    elseif i == j
        return complex(one(θ))
    else
        return complex(zero(θ))
    end
end

function R(sys::System{N}, ia, ib, θ) where N
    return Operator(sys.basis1, [Rele(ia, ib, i, j, θ) for i in 1:N, j in 1:N])
end
R2(sys::System, ia, ib, θa, θb) = R(sys, ia, ib, θa) ⊗ R(sys, ia, ib, θb)

function ratio_angle(r, ε=0)
    return 2 * asin(sqrt(r)) * (1 + ε)
end
R2_ratio(sys::System, ia, ib, r, εa, εb) =
    R2(sys, ia, ib, ratio_angle(r, εa), ratio_angle(r, εb))

apply(ρ, op) = op * ρ * op'

struct SysWrapper{Sys,Ta,Tb}
    sys::Sys
    εa::Ta
    εb::Tb
end
@inline SysWrapper(sys, εs) = SysWrapper(sys, εs[1], εs[2])
@inline R2_ratio(sys::SysWrapper, ia, ib, r) =
    R2_ratio(sys.sys, ia, ib, r, sys.εa, sys.εb)

const sys2 = System{2}()
const sys3 = System{3}()
const sys4 = System{4}()

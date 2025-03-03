#!/usr/bin/julia

using QuantumOptics
using PyPlot

using NaCsSim.Master

struct ThreeLevelData
end
const ThreeLevel = Master.SystemCoherent{Matrix{ComplexF64},3,ThreeLevelData}

struct Drive{F}
    f::F
end

function Master.init!(drive::Drive, sys::ThreeLevel)
    fill!(sys.op, 0)
end

@inline function Master.update!(drive::Drive, sys::ThreeLevel, t)
    Ω12, Ω23, Ω13 = drive.f(t)
    @inbounds begin
        op = sys.op
        op[1, 2] = Ω12 / 2
        op[2, 1] = conj(Ω12 / 2)
        op[2, 3] = Ω23 / 2
        op[3, 2] = conj(Ω23 / 2)
        op[1, 3] = Ω13 / 2
        op[3, 1] = conj(Ω13 / 2)
    end
    return
end
ThreeLevel() = Master.SystemCoherent{Matrix{ComplexF64},3}(ThreeLevelData())

function evolve(sys::ThreeLevel, dri::F, ψ0, tspan; kws...) where F
    basis = NLevelBasis(3)
    function fout(t, xs)
        n = sqrt(sum(abs2(x) for x in xs))
        return Ket(basis, [x / n for x in xs])
    end
    return Master.evolve(Drive(dri), sys, ψ0.data, tspan; fout=fout, kws...)
end

function evolve_dressed(sys::ThreeLevel, start0,
                        δdress, Ωdress, δshelve, fΩshelve::F, tspan;
                        kws...) where F
    function dri(t)
        return Ωdress * cis(δdress * t), fΩshelve(t) * cis(δshelve * t), 0
    end
    dp = δdress / sqrt(δdress^2 + Ωdress^2)
    if start0
        a1 = sqrt((1 - dp) / 2)
        a2 = sqrt((1 + dp) / 2)
    else
        a1 = sqrt((1 + dp) / 2)
        a2 = -sqrt((1 - dp) / 2)
    end
    ψinit = ComplexF64[a1, a2, 0]
    return evolve(sys, dri, Ket(NLevelBasis(3), ψinit), tspan; kws...)
end

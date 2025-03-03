#!/usr/bin/julia

using QuantumOptics
using PyPlot

using NaCsSim.Master

struct TwoLevelData
end
const TwoLevel = Master.SystemCoherent{Matrix{ComplexF64},2,TwoLevelData}

struct Drive{F}
    f::F
end

function Master.init!(drive::Drive, sys::TwoLevel)
    fill!(sys.op, 0)
end

@inline function Master.update!(drive::Drive, sys::TwoLevel, t)
    x, y, z = drive.f(t)
    @inbounds begin
        op = sys.op
        op[1, 1] = z
        op[2, 2] = -z
        op[1, 2] = complex(x, y)
        op[2, 1] = complex(x, -y)
    end
    return
end
TwoLevel() = Master.SystemCoherent{Matrix{ComplexF64},2}(TwoLevelData())
const basis_2 = SpinBasis(1//2)

function evolve(sys::TwoLevel, dri::F, ψ0, tspan; kws...) where F
    function fout(t, xs)
        n = sqrt(sum(abs2(x) for x in xs))
        return Ket(basis_2, [x / n for x in xs])
    end
    return Master.evolve(Drive(dri), sys, ψ0.data, tspan; fout=fout, kws...)
end

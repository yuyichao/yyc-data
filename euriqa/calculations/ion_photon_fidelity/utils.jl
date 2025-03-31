#!/usr/bin/julia

using NLopt
using JuMP

struct FidelityCalculator{N}
    m::Model
    obj::NonlinearExpr
    phase_vars::Vector{VariableRef}
    real_off::Vector{VariableRef}
    imag_off::Vector{VariableRef}

    function FidelityCalculator{N}() where N
        m = Model(NLopt.Optimizer)
        set_attribute(m, "algorithm", :LD_MMA)
        phase_vars = @variable(m, [1:N - 1], start=0.0)
        phases = [0; phase_vars]
        cos_phases = cos.(phases)
        sin_phases = sin.(phases)
        real_off = @variable(m, [1:(N - 1) * N ÷ 2])
        imag_off = @variable(m, [1:(N - 1) * N ÷ 2])
        k = 0
        obj = 0
        for i in 1:(N - 1)
            c1 = cos_phases[i]
            s1 = sin_phases[i]
            for j in (i + 1):N
                k += 1
                c2 = cos_phases[j]
                s2 = sin_phases[j]

                c = c1 * c2 + s1 * s2
                s = s1 * c2 - c1 * s2
                obj += real_off[k] * c - imag_off[k] * s
            end
        end
        @objective(m, Max, obj)
        return new{N}(m, obj, phase_vars, real_off, imag_off)
    end
end

function (calc::FidelityCalculator{N})(vals::Vararg{Any,N2}) where {N,N2}
    @assert N2 == (N - 1) * N
    for i in 1:N2 ÷ 2
        fix(calc.real_off[i], vals[i])
        fix(calc.imag_off[i], vals[i + N2 ÷ 2])
    end
    JuMP.optimize!(calc.m)
    return (value(calc.obj) * 2 + N) / (N * N)
end

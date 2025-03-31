#!/usr/bin/julia

using NLopt
using Ipopt
using JuMP
using FiniteDiff
using StaticArrays

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

function get_finite_grad(func)
    return function (g, args::Vararg{Any,N}) where N
        FiniteDiff.finite_difference_gradient!(g, x->func(x...),
                                               MVector{N,Float64}(args...),
                                               Val(:forward))
        return
    end
end

function complex_mul((r1, i1), (r2, i2))
    return r1 * r2 - i1 * i2, r1 * i2 + i1 * r2
end

function determinant_trivial(Mr, Mi)
    N = size(Mr, 1)
    if N == 1
        return Mr[1, 1], Mi[1, 1]
    end
    r = complex_mul((Mr[1, 1], Mi[1, 1]),
                    determinant_trivial(Mr[2:N, 2:N], Mi[2:N, 2:N]))
    for i in 2:N
        ele = complex_mul((Mr[1, i], Mi[1, i]),
                          determinant_trivial(Mr[2:N, [1:(i - 1); (i + 1):N]],
                                              Mi[2:N, [1:(i - 1); (i + 1):N]]))
        if i % 2 == 1
            r = r .+ ele
        else
            r = r .- ele
        end
    end
    return r
end

function add_pos_constraints(m, ρr, ρi)
    N = size(ρr, 1)
    if N == 1
        return
    end
    add_pos_constraints(m, ρr[1:N - 1, 1:N - 1], ρi[1:N - 1, 1:N - 1])
    @constraint(m, determinant_trivial(ρr, ρi)[1] >= 0)
end

function get_real_imag_matrix(M)
    N = size(M, 1)
    return ([(i <= j ? M[i, j] : M[j, i]) for i in 1:N, j in 1:N],
            [(i == j ? 0 : (i < j ? M[j, i] : -M[i, j])) for i in 1:N, j in 1:N])
end

struct IonIonModel{N}
    m::Model
    fcalc::FidelityCalculator{N}
    ρ1r::Matrix{VariableRef}
    ρ1i::Matrix{Any}
    ρ2r::Matrix{VariableRef}
    ρ2i::Matrix{Any}
    obj::VariableRef
    function IonIonModel{N}() where N
        fcalc = FidelityCalculator{N}()
        m = Model(Ipopt.Optimizer)
        @variable(m, ρ1[1:N, 1:N])
        @variable(m, ρ2[1:N, 1:N])
        for i in 1:N
            set_lower_bound(ρ1[i, i], 0)
            set_lower_bound(ρ2[i, i], 0)
            set_upper_bound(ρ1[i, i], 1)
            set_upper_bound(ρ2[i, i], 1)
        end
        rate1 = sum(ρ1[i, i] for i in 1:N)
        rate2 = sum(ρ2[i, i] for i in 1:N)
        ρ1r, ρ1i = get_real_imag_matrix(ρ1)
        ρ2r, ρ2i = get_real_imag_matrix(ρ2)
        add_pos_constraints(m, ρ1r, ρ1i)
        add_pos_constraints(m, ρ2r, ρ2i)

        r1 = VariableRef[]
        i1 = VariableRef[]
        r2 = VariableRef[]
        i2 = VariableRef[]
        for i in 1:(N - 1)
            for j in i + 1:N
                push!(r1, ρ1[i, j])
                push!(i1, ρ1[j, i])
                push!(r2, ρ2[i, j])
                push!(i2, ρ2[j, i])
            end
        end
        gradf = get_finite_grad(fcalc)
        @operator(m, ffunc, N * (N - 1), (x...)->fcalc(x...), gradf)
        f1 = @expression(m, ffunc(r1..., i1...) / rate1)
        f2 = @expression(m, ffunc(r2..., i2...) / rate2)
        @variable(m, f)
        @constraint(m, f >= f1)
        @constraint(m, f >= f2)
        @objective(m, Min, f)
        return new{N}(m, fcalc, ρ1r, ρ1i, ρ2r, ρ2i, f)
    end
end

function calc_mid(lb, ub)
    if lb === nothing
        return ub
    elseif ub === nothing
        return lb
    else
        return (lb + ub) / 2
    end
end

function constraint_pair!(model::IonIonModel, i, j;
                          rate_lb=nothing, rate_ub=nothing,
                          fid_lb=nothing, fid_ub=nothing)
    @assert i != j
    rate = model.ρ1r[i, i] * model.ρ2r[j, j] + model.ρ2r[i, i] * model.ρ1r[j, j]
    if rate_lb !== nothing
        @constraint(model.m, rate >= rate_lb)
    end
    if rate_ub !== nothing
        @constraint(model.m, rate <= rate_ub)
    end
    rate_mid = calc_mid(rate_lb, rate_ub)
    if rate_mid !== nothing
        set_start_value(model.ρ1r[i, i], sqrt(rate_mid / 2))
        set_start_value(model.ρ2r[i, i], sqrt(rate_mid / 2))
        set_start_value(model.ρ1r[j, j], sqrt(rate_mid / 2))
        set_start_value(model.ρ2r[j, j], sqrt(rate_mid / 2))
    end

    off_r, off_i = complex_mul((model.ρ1r[i, j], model.ρ1i[i, j]),
                               (model.ρ2r[i, j], model.ρ2i[i, j]))
    off2 = off_r^2 + off_i^2
    if fid_lb !== nothing
        @assert fid_lb >= 0.5
        @constraint(model.m, off2 >= ((fid_lb - 0.5) * rate)^2)
    end
    if fid_ub !== nothing
        @constraint(model.m, off2 <= ((fid_ub - 0.5) * rate)^2)
    end
    fid_mid = calc_mid(fid_lb, fid_ub)
    if fid_mid !== nothing && rate_mid !== nothing
        set_start_value(model.ρ1r[i, j], ((fid_mid - 0.5) * rate_mid))
        set_start_value(model.ρ1i[i, j], 0)
        set_start_value(model.ρ2r[i, j], ((fid_mid - 0.5) * rate_mid))
        set_start_value(model.ρ2i[i, j], 0)
    end
end

function min_fidelity!(model::IonIonModel)
    JuMP.optimize!(model.m)
    return value(calc.obj)
end

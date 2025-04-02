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
        imag_off = @variable(m, [1:(N - 1) * N ÷ 2 - N + 1])
        obj = real_off[1] * cos_phases[2]
        k = 1
        for j in 3:N
            k += 1
            obj += real_off[k] * cos_phases[j]
        end
        for i in 2:(N - 1)
            c1 = cos_phases[i]
            s1 = sin_phases[i]
            for j in (i + 1):N
                k += 1
                c2 = cos_phases[j]
                s2 = sin_phases[j]

                c = c1 * c2 + s1 * s2
                s = s1 * c2 - c1 * s2
                obj += real_off[k] * c - imag_off[k - N + 1] * s
            end
        end
        @objective(m, Max, obj)
        return new{N}(m, obj, phase_vars, real_off, imag_off)
    end
end

function (calc::FidelityCalculator{N})(vals::Vararg{Any,N2}) where {N,N2}
    @assert N2 == (N - 2) * N + 1
    N_ = (N - 1) * N ÷ 2
    for i in 1:N_
        fix(calc.real_off[i], vals[i])
        if i > N - 1
            fix(calc.imag_off[i - N + 1], vals[i - N + 1 + N_])
        end
    end
    JuMP.optimize!(calc.m)
    return value(calc.obj)
end

function get_finite_grad(func)
    return function (g, args::Vararg{Any,N}) where N
        FiniteDiff.finite_difference_gradient!(g, x->func(x...),
                                               MVector{N,Float64}(args...),
                                               Val(:forward))
        return
    end
end

function rmul(r1, r2)
    if iszero(r1) || iszero(r2)
        return 0
    elseif isa(r1, Number) && isone(r1)
        return r2
    elseif isa(r1, Number) && isone(r2)
        return r1
    end
    return r1 * r2
end

function radd(r1, r2)
    if iszero(r1)
        return r2
    elseif iszero(r2)
        return r1
    end
    return r1 + r2
end

function rsub(r1, r2)
    if iszero(r1)
        return -r2
    elseif iszero(r2)
        return r1
    end
    return r1 - r2
end

function cmul((r1, i1), (r2, i2))
    return (rsub(rmul(r1, r2), rmul(i1, i2)),
            radd(rmul(r1, i2), rmul(i1, r2)))
end

function cadd((r1, i1), (r2, i2))
    return radd(r1, r2), radd(i1, i2)
end

function csub((r1, i1), (r2, i2))
    return rsub(r1, r2), rsub(i1, i2)
end

function determinant_trivial(Mr, Mi)
    N = size(Mr, 1)
    if N == 1
        return Mr[1, 1], Mi[1, 1]
    end
    r = cmul((Mr[1, 1], Mi[1, 1]),
             determinant_trivial(Mr[2:N, 2:N], Mi[2:N, 2:N]))
    for i in 2:N
        ele = cmul((Mr[1, i], Mi[1, i]),
                   determinant_trivial(Mr[2:N, [1:(i - 1); (i + 1):N]],
                                       Mi[2:N, [1:(i - 1); (i + 1):N]]))
        if i % 2 == 1
            r = cadd(r, ele)
        else
            r = csub(r, ele)
        end
    end
    return r
end

function add_pos_constraints(m, constraints, ρr, ρi)
    N = size(ρr, 1)
    for i in 1:N - 1
        for j in i + 1:N
            ex = determinant_trivial(ρr[[i, j], [i, j]], ρi[[i, j], [i, j]])[1]
            push!(constraints, ex)
            @NLconstraint(m, ex >= 0)
        end
    end
    for i in 1:N - 2
        for j in i + 1:N - 1
            for k in j + 1:N
                ex = determinant_trivial(ρr[[i, j, k], [i, j, k]],
                                         ρi[[i, j, k], [i, j, k]])[1]
                push!(constraints, ex)
                @NLconstraint(m, ex >= 0)
            end
        end
    end
    for i in 4:N
        ex = determinant_trivial(ρr[1:i, 1:i], ρi[1:i, 1:i])[1]
        push!(constraints, ex)
        @NLconstraint(m, ex >= 0)
    end
end

function create_density_matrix(m, name, N)
    ρr = Matrix{Any}(undef, N, N)
    ρi = Matrix{Any}(undef, N, N)
    fid_args_r = []
    fid_args_i = []
    for i in 1:N
        for j in i:N
            er = @variable(m, base_name="$(name)r[$i, $j]", start = 0)
            ρr[i, j] = er
            set_upper_bound(er, 1)
            if i == j
                set_lower_bound(er, 0)
                ρi[i, j] = 0
                continue
            end
            set_lower_bound(er, -1)
            push!(fid_args_r, er)
            ρr[j, i] = er
            if i == 1
                ρi[i, j] = 0
                ρi[j, i] = 0
                continue
            end
            ei = @variable(m, base_name="$(name)i[$i, $j]", start = 0)
            set_lower_bound(ei, -1)
            set_upper_bound(ei, 1)
            ρi[i, j] = ei
            ρi[j, i] = -ei
            push!(fid_args_i, ei)
        end
    end
    return ρr, ρi, [fid_args_r; fid_args_i]
end

struct IonIonModel{N}
    m::Model
    fcalc::FidelityCalculator{N}
    constraints::Vector{Any}
    ρ1r::Matrix{Any}
    ρ1i::Matrix{Any}
    ρ2r::Matrix{Any}
    ρ2i::Matrix{Any}
    obj::Any
    function IonIonModel{N}(m::Model) where N
        fcalc = FidelityCalculator{N}()
        ρ1r, ρ1i, fid_args1 = create_density_matrix(m, "ρ1", N)
        ρ2r, ρ2i, fid_args2 = create_density_matrix(m, "ρ2", N)
        rate1 = @NLexpression(m, sum(ρ1r[i, i] for i in 1:N))
        rate2 = @NLexpression(m, sum(ρ2r[i, i] for i in 1:N))
        constraints = []
        add_pos_constraints(m, constraints, ρ1r, ρ1i)
        add_pos_constraints(m, constraints, ρ2r, ρ2i)
        gradf = get_finite_grad(fcalc)
        register(m, :ffunc, N * (N - 2) + 1, (x...)->fcalc(x...), gradf, autodiff=false)
        # @operator(m, ffunc, N * (N - 2) + 1, (x...)->fcalc(x...), gradf)
        f1 = @NLexpression(m, (ffunc(fid_args1...) / rate1 * 2 + 1) / N)
        f2 = @NLexpression(m, (ffunc(fid_args2...) / rate2 * 2 + 1) / N)
        f = @NLexpression(m, max(f1, f2))
        @NLobjective(m, Min, f)
        return new{N}(m, fcalc, constraints, ρ1r, ρ1i, ρ2r, ρ2i, f)
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

function rate_expr(model::IonIonModel, i, j)
    return model.ρ1r[i, i] * model.ρ2r[j, j] + model.ρ2r[i, i] * model.ρ1r[j, j]
end
function fid_expr(model::IonIonModel, i, j)
    off_r, off_i = cmul((model.ρ1r[i, j], model.ρ1i[i, j]),
                        (model.ρ2r[i, j], model.ρ2i[i, j]))
    off2 = radd(off_r^2, off_i^2)
    return 0.5 + sqrt(off2) / rate_expr(model, i, j)
end

function constraint_pair!(model::IonIonModel, i, j;
                          rate_lb=nothing, rate_ub=nothing,
                          fid_lb=nothing, fid_ub=nothing)
    @assert i != j
    if i > j
        i, j = j, i
    end
    rate = rate_expr(model, i, j)
    if rate_lb !== nothing
        @NLconstraint(model.m, rate >= rate_lb)
    end
    if rate_ub !== nothing
        @NLconstraint(model.m, rate <= rate_ub)
    end
    rate_mid = calc_mid(rate_lb, rate_ub)
    if rate_mid !== nothing
        set_start_value(model.ρ1r[i, i], sqrt(rate_mid / 2))
        set_start_value(model.ρ2r[i, i], sqrt(rate_mid / 2))
        set_start_value(model.ρ1r[j, j], sqrt(rate_mid / 2))
        set_start_value(model.ρ2r[j, j], sqrt(rate_mid / 2))
    end

    off_r, off_i = cmul((model.ρ1r[i, j], model.ρ1i[i, j]),
                        (model.ρ2r[i, j], model.ρ2i[i, j]))
    off2 = radd(off_r^2, off_i^2)
    if fid_lb !== nothing
        @assert fid_lb >= 0.5
        @NLconstraint(model.m, off2 >= ((fid_lb - 0.5) * rate)^2)
    end
    if fid_ub !== nothing
        @NLconstraint(model.m, off2 <= ((fid_ub - 0.5) * rate)^2)
    end
    fid_mid = calc_mid(fid_lb, fid_ub)
    if fid_mid !== nothing && rate_mid !== nothing
        set_start_value(model.ρ1r[i, j], sqrt((fid_mid - 0.5) * rate_mid))
        set_start_value(model.ρ2r[i, j], sqrt((fid_mid - 0.5) * rate_mid))
    end
end

function min_fidelity!(model::IonIonModel)
    JuMP.optimize!(model.m)
    return value(model.obj)
end

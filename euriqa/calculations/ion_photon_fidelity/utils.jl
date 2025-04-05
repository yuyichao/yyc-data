#!/usr/bin/julia

using NLopt
using Ipopt
using JuMP
using FiniteDiff
using StaticArrays
using LinearAlgebra

struct FidelityCalculator{N}
    m::Model
    obj::NonlinearExpr
    phase_vars::Vector{VariableRef}
    real_off::Vector{VariableRef}
    imag_off::Vector{VariableRef}

    function FidelityCalculator{N}() where N
        m = Model(NLopt.Optimizer)
        set_attribute(m, "algorithm", :LD_TNEWTON_PRECOND_RESTART)
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

function add_pos_constraints(m, ρr, ρi)
    N = size(ρr, 1)
    for i in 1:N - 1
        for j in i + 1:N
            ex = determinant_trivial(ρr[[i, j], [i, j]], ρi[[i, j], [i, j]])[1]
            @NLconstraint(m, ex >= 0)
        end
    end
    for i in 1:N - 2
        for j in i + 1:N - 1
            for k in j + 1:N
                ex = determinant_trivial(ρr[[i, j, k], [i, j, k]],
                                         ρi[[i, j, k], [i, j, k]])[1]
                @NLconstraint(m, ex >= 0)
            end
        end
    end
    for i in 4:N
        ex = determinant_trivial(ρr[1:i, 1:i], ρi[1:i, 1:i])[1]
        @NLconstraint(m, ex >= 0)
    end
end

function create_density_matrix(m, name, N)
    ρr = Matrix{Any}(undef, N, N)
    ρi = Matrix{Any}(undef, N, N)
    fid_args_d = []
    fid_args_r = []
    fid_args_i = []
    for i in 1:N
        for j in i:N
            er = @variable(m, base_name="$(name)r[$i, $j]", start = 0)
            ρr[i, j] = er
            set_upper_bound(er, 1)
            if i == j
                push!(fid_args_d, er)
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
    return ρr, ρi, [fid_args_d; fid_args_r; fid_args_i]
end

struct IonIonModel{N}
    m::Model
    fcalc::FidelityCalculator{N}
    ρ1r::Matrix{Any}
    ρ1i::Matrix{Any}
    ρ2r::Matrix{Any}
    ρ2i::Matrix{Any}
    obj::Any
    vars::Vector{VariableRef}
    function IonIonModel{N}(m::Model) where N
        fcalc = FidelityCalculator{N}()
        ρ1r, ρ1i, fid_args1 = create_density_matrix(m, "ρ1", N)
        ρ2r, ρ2i, fid_args2 = create_density_matrix(m, "ρ2", N)
        rate1 = @NLexpression(m, sum(ρ1r[i, i] for i in 1:N))
        rate2 = @NLexpression(m, sum(ρ2r[i, i] for i in 1:N))
        add_pos_constraints(m, ρ1r, ρ1i)
        add_pos_constraints(m, ρ2r, ρ2i)
        function full_fid(x...)
            diag1 = x[1:N]
            off1 = x[N + 1:N * (N - 1) + 1]
            diag2 = x[N * (N - 1) + 2:N * N + 1]
            off2 = x[N * N + 2:end]
            return (max(fcalc(off1...) / sum(diag1),
                        fcalc(off2...) / sum(diag2)) * 2 + 1) / N
        end
        gradf = get_finite_grad(full_fid)
        register(m, :ffunc, (N * (N - 1) + 1) * 2, full_fid, gradf, autodiff=false)
        f = @NLexpression(m, ffunc(fid_args1..., fid_args2...))
        @NLobjective(m, Min, f)
        return new{N}(m, fcalc, ρ1r, ρ1i, ρ2r, ρ2i, f,
                      [fid_args1..., fid_args2...])
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

function constraint_pair!(model::IonIonModel, i, j; rate_lb, rate_ub, fid_lb, fid_ub)
    @assert i != j
    if i > j
        i, j = j, i
    end
    rate = rate_expr(model, i, j)
    @NLconstraint(model.m, rate >= rate_lb)
    @NLconstraint(model.m, rate <= rate_ub)
    rate_mid = (rate_lb + rate_ub) / 2
    set_start_value(model.ρ1r[i, i], sqrt(rate_mid / 2))
    set_start_value(model.ρ2r[i, i], sqrt(rate_mid / 2))
    set_start_value(model.ρ1r[j, j], sqrt(rate_mid / 2))
    set_start_value(model.ρ2r[j, j], sqrt(rate_mid / 2))

    off_r, off_i = cmul((model.ρ1r[i, j], model.ρ1i[i, j]),
                        (model.ρ2r[i, j], model.ρ2i[i, j]))
    off2 = radd(off_r^2, off_i^2)
    @assert fid_lb >= 0.5
    @NLconstraint(model.m, off2 >= ((fid_lb - 0.5) * rate)^2)
    @NLconstraint(model.m, off2 <= ((fid_ub - 0.5) * rate)^2)
    set_start_value(model.ρ1r[i, j], sqrt((fid_lb - 0.5) * rate_mid))
    set_start_value(model.ρ2r[i, j], sqrt((fid_lb - 0.5) * rate_mid))
end

function min_fidelity!(model::IonIonModel)
    JuMP.optimize!(model.m)
    return value(model.obj)
end

struct IonIonConstraints{N}
    bounds::Dict{NTuple{2,Int},NTuple{4,Float64}}
    function IonIonConstraints{N}() where N
        return new{N}(Dict{NTuple{2,Int},NTuple{4,Float64}}())
    end
end
function constraint_pair!(c::IonIonConstraints, i, j; rate_lb, rate_ub, fid_lb, fid_ub)
    @assert i != j
    if i > j
        i, j = j, i
    end
    c.bounds[(i, j)] = (rate_lb, rate_ub, fid_lb, fid_ub)
    return
end

function rand_find_min_fidelity(fcalc, vals, ::Val{N}) where N
    diag_offset = 0
    offr_offset = N
    offr2_offset = 2N - 1
    offi2_offset = N * (N + 1) ÷ 2
    total_vals = N * (N - 1) + 1
    @assert total_vals == length(vals)
    M = MMatrix{N,N,ComplexF64}(undef)
    fid_args_r1 = MVector{N - 1,Float64}(undef)
    for i in 1:N
        M[i, i] = vals[i]
        if i != 1
            M[1, i] = M[i, 1] = fid_args_r1[i - 1] = vals[i - 1 + offr_offset]
        end
    end
    fid_args_a2 = MVector{(N - 1) * (N - 2) ÷ 2,Float64}(undef)
    for i in 1:length(fid_args_a2)
        fid_args_a2[i] = hypot(vals[i + offr2_offset], vals[i + offi2_offset])
    end
    fid_args_r2 = MVector{(N - 1) * (N - 2) ÷ 2,Float64}(undef)
    fid_args_i2 = MVector{(N - 1) * (N - 2) ÷ 2,Float64}(undef)
    function rand_off!()
        for i in 1:length(fid_args_a2)
            im_scale = rand() * 0.2 - 0.1
            re_scale = sqrt(1 - im_scale^2)
            fid_args_r2[i] = fid_args_a2[i] * re_scale
            fid_args_i2[i] = fid_args_a2[i] * im_scale
        end
    end
    function compute_fid()
        return fcalc(fid_args_r1..., fid_args_r2..., fid_args_i2...)
    end
    function check_posdef()
        off_idx = 0
        for i in 2:N - 1
            for j in i + 1:N
                off_idx += 1
                vr = fid_args_r2[off_idx]
                vi = fid_args_i2[off_idx]
                M[i, j] = complex(vr, vi)
                M[j, i] = complex(vr, -vi)
            end
        end
        return isposdef(M)
    end
    fid_args_r2 .= fid_args_a2
    fid_args_i2 .= 0

    opt_args = (fid_args_r1..., fid_args_r2..., fid_args_i2...)
    opt_fid = compute_fid()
    for _ in 1:20000
        rand_off!()
        if !check_posdef()
            continue
        end
        fid = compute_fid()
        if fid < opt_fid
            opt_args = (fid_args_r1..., fid_args_r2..., fid_args_i2...)
            opt_fid = fid
        end
    end
    return (opt_fid / sum(@view(vals[1:N])) * 2 + 1) / N, opt_args
end

function min_fidelity!(c::IonIonConstraints{N}) where N
    m = Model(NLopt.Optimizer)
    set_attribute(m, "algorithm", :LD_SLSQP)
    model = IonIonModel{N}(m)
    for i in 1:N - 1
        for j in (i + 1):N
            rate_lb, rate_ub, fid_lb, fid_ub =
                get(c.bounds, (i, j), (0.0, 1 / N^2, 0.0, 1.0))
            constraint_pair!(model, i, j, rate_lb=rate_lb, rate_ub=rate_ub,
                             fid_lb=fid_lb, fid_ub=fid_ub)
        end
    end
    fid1 = @time min_fidelity!(model)
    println("Initial result: $fid1")
    fcalc = FidelityCalculator{N}()
    var_vals = value.(model.vars)
    Nsingle = N * (N - 1) + 1
    @assert length(var_vals) == Nsingle * 2
    opt_fid1, opt_args1 = @time rand_find_min_fidelity(fcalc, var_vals[1:Nsingle], Val(N))
    opt_fid2, opt_args2 = @time rand_find_min_fidelity(fcalc, var_vals[Nsingle + 1:end], Val(N))
    println("Randomized result: $opt_fid1, $opt_fid2")
    set_start_value.(model.vars, var_vals)
    for i in 1:length(opt_args1)
        set_start_value(model.vars[N + i], opt_args1[i])
        set_start_value(model.vars[N + i + Nsingle], opt_args2[i])
    end
    fid2 = @time(min_fidelity!(model))
    println("Final result: $fid2")
    return fid2, model
end

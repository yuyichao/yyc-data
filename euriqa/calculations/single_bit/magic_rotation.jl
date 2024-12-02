#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using PyPlot
using NaCsPlot
using JuMP
using NLopt
using Ipopt

include("pauli_utils.jl")

function target_y(Ω, cutoff)
    return Ω > cutoff ? 0 : 1
end

function single_rotation(α, X, Z, ϕ)
    X = α * X
    A = sqrt(X^2 + Z^2)
    sincA = sinc(A / π)
    sϕ, cϕ = sincos(ϕ)
    return (cos(A), im * sincA * X * cϕ, im * sincA * X * sϕ, im * sincA * Z)
end

function multi_rotation(α, param0, params...)
    v = single_rotation(α, param0...)
    for param in params
        v = multiply_pauli(v, single_rotation(α, param...))
    end
    return v
end

function multi_rotation_xy(α, params...)
    v = multi_rotation(α, params...)
    vxy, vz = decompose_xy_z(v)
    return vxy[2], vxy[3]
end

function objective_function(α, Ω, cutoff, params...)
    x, y = multi_rotation_xy(α, params...)
    ytgt = target_y(α * Ω, cutoff)
    ydiff = imag(y) - target_y(α * Ω, cutoff)
    if abs(ydiff) > 1
        ydiff = 2 - abs(ydiff)
    end
    return abs2(x) + abs2(ydiff)
end

function repack_param3(params::Vararg{Any,N}) where N
    let N_3 = N ÷ 3
        return ntuple(i->(params[i], params[i + N_3], params[i + N_3 * 2]), N_3)
    end
end

function opt_n(model, ::Val{n}, Ω, cutoff, Xinit=pi/8, Zinit=0.1, ϕinit=pi/2;
               minimize_angle=false, allow_z=false, robust=false) where n
    function obj_func(p1...)
        params = repack_param3(p1...)
        if robust
            return (objective_function(0.98, Ω, cutoff, params...) * 1 +
                objective_function(1.0, Ω, cutoff, params...) * 5 +
                objective_function(1.02, Ω, cutoff, params...) * 1 +
                objective_function(√(2) - 0.02, Ω, cutoff, params...) * 1 +
                objective_function(√(2), Ω, cutoff, params...) * 5 +
                objective_function(√(2) + 0.02, Ω, cutoff, params...) * 1)
        else
            return (objective_function(1.0, Ω, cutoff, params...) * 5 +
                objective_function(√(2), Ω, cutoff, params...) * 5)
        end
    end

    Xv = Vector{Float64}(undef, n)
    Zv = Vector{Float64}(undef, n)
    ϕv = Vector{Float64}(undef, n)
    if model !== nothing
        empty!(model)
        X = ntuple(_->@variable(model, lower_bound=0, upper_bound=2π), n)
        set_start_value.(X, Xinit)
        if allow_z
            Z = ntuple(_->@variable(model, lower_bound=-4π, upper_bound=4π), n)
            set_start_value.(Z, Zinit)
        else
            Z = ntuple(_->0, n)
        end
        ϕ = ntuple(_->@variable(model, lower_bound=0, upper_bound=3π), n)
        set_start_value.(ϕ, ϕinit)
        P = tuple(X..., Z..., ϕ...)

        # @operator(model, obj_f, 3 * n, obj_func)
        # @objective(model, Min, obj_f(P...))

        register(model, :obj_f, 3 * n, obj_func, autodiff=true)
        if minimize_angle
            @NLobjective(model, Min, obj_f(P...) * (+(X...) + 1))
        else
            @NLobjective(model, Min, obj_f(P...))
        end
        JuMP.optimize!(model)
        Xv .= value.(X)
        if allow_z
            Zv .= value.(Z)
        else
            Zv .= 0
        end
        ϕv .= value.(ϕ)
    else
        Xv .= Xinit
        if allow_z
            Zv .= Zinit
        else
            Zv .= 0
        end
        ϕv .= ϕinit
    end
    return Xv, Zv, ϕv, obj_func(Xv..., Zv..., ϕv...)
end

function plot_xy(αs, Xs, Zs, ϕs)
    xs = similar(αs, Float64)
    ys = similar(αs, Float64)
    for (i, α) in enumerate(αs)
        x, y = multi_rotation_xy(α, zip(Xs, Zs, ϕs)...)
        xs[i] = imag(x)
        ys[i] = imag(y)
    end
    return xs, ys
end

# [1.573808755032251, 4.720886497147208, 1.5734306439589152], [0.0, 0.0, 0.0], [3.4330944508917947, 5.114712612318794, 0.5106233267300696], 0.0001138493952966524

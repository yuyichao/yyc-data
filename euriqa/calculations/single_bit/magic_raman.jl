#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using PyPlot
using NaCsPlot
using JuMP
using NLopt
using Ipopt

include("pauli_utils.jl")

function single_step(α, Ω0, δ0, rzx, t1, t2)
    Ω = α * Ω0
    δ = Ω * rzx + δ0
    X1 = Ω * t1
    Z1 = δ * t1
    Z2 = δ0 * t2
    A1 = sqrt(X1^2 + Z1^2)
    sincA1 = sinc(A1 / π)
    σ1 = (cos(A1), im * sincA1 * X1, 0, im * sincA1 * Z1)
    s2, c2 = sincos(Z2 / 2)
    σ2 = (c2, 0, 0, im * s2)
    return multiply_pauli(σ2, σ1)
end

function multi_step(α, Ω0, δ0, rzx, param0, params...)
    v = single_step(α, Ω0, δ0, rzx, param0...)
    for param in params
        v = multiply_pauli(v, single_step(α, Ω0, δ0, rzx, param...))
    end
    return v
end

function multi_step_xyz(args...)
    vxy, vz, c = decompose_xy_z(multi_step(args...))
    return vxy[2], vxy[3], vz[4], c
end

function objective_function(x0, args...)
    x, y, z, c = multi_step_xyz(args...)
    return abs2(imag(x) - x0) + abs2(y) + abs2(z)
end

function repack_param2(params::Vararg{Any,N}) where N
    let N_2 = N ÷ 2
        return ntuple(i->(params[i], params[i + N_2]), N_2)
    end
end

function opt_n(model, ::Val{n}, x0, δ0, rzx, Ω0init=-δ0 / rzx, t1init=0.3/Ω0, t2init=0.3/δ0;
               minimize_time=false) where n
    function obj_func(Ω0, p1...)
        params = repack_param2(p1...)
        return (objective_function(x0, 1 - 0.06, Ω0, δ0, rzx, params...) * 1 +
            objective_function(x0, 1 - 0.03, Ω0, δ0, rzx, params...) * 3 +
            objective_function(x0, 1, Ω0, δ0, rzx, params...) * 5 +
            objective_function(x0, 1 + 0.03, Ω0, δ0, rzx, params...) * 3 +
            objective_function(x0, 1 + 0.06, Ω0, δ0, rzx, params...) * 1)
    end

    t1v = Vector{Float64}(undef, n)
    t2v = Vector{Float64}(undef, n)
    if model !== nothing
        empty!(model)
        Ω0 = @variable(model, lower_bound=0)
        set_start_value.(Ω0, Ω0init)
        t1 = ntuple(_->@variable(model, lower_bound=0, upper_bound=6π / abs(δ0)), n)
        set_start_value.(t1, t1init)
        t2 = ntuple(_->@variable(model, lower_bound=0, upper_bound=4π / abs(δ0)), n)
        set_start_value.(t2, t2init)
        P = tuple(t1..., t2...)

        # @operator(model, obj_f, 3 * n, obj_func)
        # @objective(model, Min, obj_f(P...))

        register(model, :obj_f, 2 * n + 1, obj_func, autodiff=true)
        if minimize_time
            @NLobjective(model, Min, obj_f(Ω0, P...) * (+(t1..., t2...) + 1 / δ0))
        else
            @NLobjective(model, Min, obj_f(Ω0, P...))
        end
        JuMP.optimize!(model)
        Ω0v = value(Ω0)
        t1v .= value.(t1)
        t2v .= value.(t2)
    else
        Ω0v = Ω0init
        t1v .= t1init
        t2v .= t2init
    end
    return Ω0v, t1v, t2v, obj_func(Ω0v, t1v..., t2v...)
end

function plot_xyz(αs, Ω0, δ0, rzx, t1s, t2s)
    xs = similar(αs, Float64)
    ys = similar(αs, Float64)
    zs = similar(αs, Float64)
    cs = similar(αs, Float64)
    for (i, α) in enumerate(αs)
        x, y, z, c = multi_step_xyz(α, Ω0, δ0, rzx, zip(t1s, t2s)...)
        xs[i] = imag(x)
        ys[i] = imag(y)
        zs[i] = imag(z)
        cs[i] = real(c)
    end
    return xs, ys, zs, cs
end

# For pi/2, δ0=-130, rzx=13/70
# 281.0215486884919, [0.004092442979960276, 0.007176534729743897, 0.005500575584942863], [0.03197704573604141, 0.014138893848454709, 0.021129655940143557], 1.02e-5

# For pi, δ0=-130, rzx=13/70
# 605.2517681505016, [0.18261584293714625, 0.16932841675170376, 0.1693849604657854], [0.04771810080803878, 0.0003627750891943768, 0.0009106648968329537], 0.000126

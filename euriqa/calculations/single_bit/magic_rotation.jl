#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using PyPlot
using NaCsPlot
using JuMP
using NLopt
using Ipopt

include("pauli_utils.jl")

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

function multi_rotation_xyz(α, params...)
    v = multi_rotation(α, params...)
    vxy, vz, c = decompose_xy_z(v)
    return vxy[2], vxy[3], vz[4], c
end

function objective_function(α, Ω, cutoff, inverse, params...)
    x, y, z, c = multi_rotation_xyz(α, params...)
    if (α * Ω > cutoff) ⊻ inverse
        return abs2(x) + abs2(y) + abs2(z) + abs2(c + 1)
    else
        return abs2(x) + abs2(abs(y) - 1)
    end
end

function repack_param3(params::Vararg{Any,N}) where N
    let N_3 = N ÷ 3
        return ntuple(i->(params[i], params[i + N_3], params[i + N_3 * 2]), N_3)
    end
end

function opt_n(model, ::Val{n}, Ω, cutoff, Xinit=pi/8, Zinit=0.1, ϕinit=pi/2;
               minimize_angle=false, allow_z=false, robust=false,
               inverse=false) where n
    function obj_func(p1...)
        params = repack_param3(p1...)
        v0 = 1.0
        v1 = cutoff * 2 - 1
        if robust
            return (objective_function(v0 - 0.02, Ω, cutoff, inverse, params...) * 1 +
                objective_function(v0, Ω, cutoff, inverse, params...) * 5 +
                objective_function(v0 + 0.02, Ω, cutoff, inverse, params...) * 1 +
                objective_function(v1 - 0.02, Ω, cutoff, inverse, params...) * 1 +
                objective_function(v1, Ω, cutoff, inverse, params...) * 5 +
                objective_function(v1 + 0.02, Ω, cutoff, inverse, params...) * 1)
        else
            return (objective_function(v0, Ω, cutoff, inverse, params...) * 5 +
                objective_function(v1, Ω, cutoff, inverse, params...) * 5)
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

function plot_xyz(αs, Xs, Zs, ϕs)
    xs = similar(αs, Float64)
    ys = similar(αs, Float64)
    zs = similar(αs, Float64)
    cs = similar(αs, Float64)
    for (i, α) in enumerate(αs)
        x, y, z, c = multi_rotation_xyz(α, zip(Xs, Zs, ϕs)...)
        xs[i] = imag(x)
        ys[i] = imag(y)
        zs[i] = imag(z)
        cs[i] = real(c)
    end
    return xs, ys, zs, cs
end

# [1.573808755032251, 4.720886497147208, 1.5734306439589152], [0.0, 0.0, 0.0], [3.4330944508917947, 5.114712612318794, 0.5106233267300696], 0.0001138493952966524
# [2.2256258753645315, 2.2238315316814696, 4.4452007353124365], [0.0, 0.0, 0.0], [4.0718830352918784, 6.499580993261763, 7.954649435204307]
# [4.150729797709635, 3.7895211205522674, 4.338291742705923], [-1.807628358753955, 3.07132682114769, 0.036334424920617926], [2.686170721467204, 0.9383493612406526, 5.123962927102851], 0.0006416825788214378
# [3.1392191703957866, 6.281845773039094, 6.283185367330217], [0.0, 0.0, 0.0], [3.8391734302197564, 1.9553083999084317, 5.601573269655488], 0.00034835886361337515
# [2.0907410201241396, 4.377400371525149, 4.291498651282052], [0.0, 0.0, 0.0], [3.5187705950221755, 5.00101379482833, 7.115251542515617], 0.012262573226535699
# [4.316028274228701, 6.110790073787895, 4.058180619416489], [-1.332004120147513, 3.6381574305945814, -2.4772392440773476], [2.7420115615339076, 5.132624320244106, 1.0755740883569926], 9.017731555765204e-5

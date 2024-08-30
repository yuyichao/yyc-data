#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using PyPlot
using NaCsPlot
using JuMP
using NLopt
using Ipopt

include("pauli_utils.jl")

function target_y(Ω, δ)
    return sqrt((1 - δ / sqrt(δ^2 + Ω^2)) / 2)
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

function objective_function(α, Ω, δ, params...)
    x, y = multi_rotation_xy(α, params...)
    return abs2(x) + abs2(imag(y) - target_y(α * Ω, δ))
end

function repack_param3(params::Vararg{Any,N}) where N
    let N_3 = N ÷ 3
        return ntuple(i->(params[i], params[i + N_3], params[i + N_3 * 2]), N_3)
    end
end

function opt_n(model, ::Val{n}, Ω, δ, Xinit=pi/8, Zinit=0.1, ϕinit=pi/2;
               minimize_angle=false) where n
    function obj_func(p1...)
        params = repack_param3(p1...)
        return (objective_function(0.005, Ω, δ, params...) * 10 +
            objective_function(0.01, Ω, δ, params...) * 10 +
            objective_function(0.02, Ω, δ, params...) * 1 +
            objective_function(0.9, Ω, δ, params...) * 3 +
            objective_function(0.95, Ω, δ, params...) * 3 +
            objective_function(0.975, Ω, δ, params...) * 3 +
            objective_function(1.0, Ω, δ, params...) * 3 +
            objective_function(1.025, Ω, δ, params...) * 2 +
            objective_function(1.05, Ω, δ, params...) * 1)
    end

    Xv = Vector{Float64}(undef, n)
    Zv = Vector{Float64}(undef, n)
    ϕv = Vector{Float64}(undef, n)
    if model !== nothing
        empty!(model)
        X = ntuple(_->@variable(model, lower_bound=0, upper_bound=2π), n)
        set_start_value.(X, Xinit)
        Z = ntuple(_->@variable(model, lower_bound=-4π, upper_bound=4π), n)
        set_start_value.(Z, Zinit)
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
        Zv .= value.(Z)
        ϕv .= value.(ϕ)
    else
        Xv .= Xinit
        Zv .= Zinit
        ϕv .= ϕinit
    end
    return Xv, Zv, ϕv, obj_func(Xv..., Zv..., ϕv...)
end

# sum(X) = 9.291310413820183
# objective = 1.3847397116613161e-6
# [4.102159441616568, 2.199318322737191, 2.989832649466424]
# [2.031387365486914, 0.13492361591220967, -0.00899924603161387]
# [2.524307369373361, 5.665893863789733, 3.940672111232241]

# Plot 1
# sum(X) = 9.272810270360356
# objective = 1.3859881161318411e-6
# [4.093345622593074, 2.1909427292438184, 2.9885219185234635]
# [2.0268537539777705, 0.13481328004661897, -0.009067422386459207]
# [2.5324832056395623, 5.6657203351806125, 3.940979719209872]

# sum(X) = 13.006179134451575
# objective = 8.102088335449348e-7
# [6.279694773817025, 1.9183038928921166, 2.649912355364408, 2.1582681123780243],
# [1.6380494351314936, -0.34877001533828, 0.0852699937152111, -1.327486449705353],
# [4.820205472294657, 1.6974691988479882, 1.6303756977163717, 0.9212896440280046]

# sum(Xv) = 8.473077303595474
# objective = 1.6110420903155428e-6
# [3.0511873609811815, 1.6863317648627598, 3.735558177751533]
# [0.20177893195957936, 2.931257488639002, -1.1391876473354867]
# [4.69721890717192, 1.9193738094854447e-6, 2.106245748078913]

# sum(Xv) = 8.471741326340982
# objective = 1.6111512584032458e-6
# [3.0506445592948133, 1.6861668598475104, 3.734929907198659]
# [0.20148282157496725, 2.93121133991079, -1.1390151354518243]
# [0.7475697026895737, 0, 0.3352218813583657] .* 2pi

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


# X[3.0506445592948133, 1.6861668598475104, 3.734929907198659]
# Z[0.20148282157496725, 2.93121133991079, -1.1390151354518243]
# axis: [0.7475697026895737, 0, 0.3352218813583657]
# δ/Ω: [0.0660459839416828, 1.7383874690645242, -0.3049629213272517]
# angle: [0.9710503224563323, 0.5367235812449408, 1.1888651136648412]

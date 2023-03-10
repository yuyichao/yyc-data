#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using PyPlot
using NaCsPlot
using JuMP
using NLopt

include("pauli_utils.jl")

function equivalent_xy_angle(M)
    σ_xy, = decompose_xy_z(M)
    return decompose_xy_rot(σ_xy)[1]
end

function gen_two_rotations(a2, c2, x_a)
    a = a2 / 2
    c = c2 / 2
    b = sqrt(a^2 - c^2)

    y = b * sqrt(1 - x_a^2)
    x = x_a * a

    x1 = x + c
    y1 = y
    x2 = c - x
    y2 = -y
    return (xy_to_polar(x1, y1), xy_to_polar(x2, y2))
end

function gen_two_rotations_M(a2, c2, x_a)
    r1, r2 = gen_two_rotations(a2, c2, x_a)
    return (compose_xy_rot(r1...), compose_xy_rot(r2...))
end

function gen_target_2π(a2, c2)
    return function (x)
        M1, M2 = gen_two_rotations_M(a2, c2, x)
        return equivalent_xy_angle(M1 * M2)
    end
end

function gen_target_π(a2, c2)
    return function (x)
        (ψ1, θ1), (ψ2, θ2) = gen_two_rotations(a2, c2, x)
        M1 = compose_xy_rot(ψ1, θ1)
        M2 = compose_xy_rot(ψ2, -θ2)
        M = M1 * compose_sigmas(0, 1, 0, 0) * M2
        return equivalent_xy_angle(M)
    end
end

function opt_angle_2π(a2, c2)
    f = gen_target_2π(a2, c2)
    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    register(model, :f, 1, f, autodiff=true)
    @variable(model, 0 <= x <= 1)
    @NLobjective(model, Max, f(x))
    JuMP.optimize!(model)
    xv = value(x)
    return xv, f(xv)
end

function opt_angle_π(a2, c2)
    f = gen_target_π(a2, c2)
    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    register(model, :f, 1, f, autodiff=true)
    @variable(model, 0 <= x <= 1)
    @NLobjective(model, Max, f(x))
    JuMP.optimize!(model)
    xv = value(x)
    return xv, f(xv)
end

function opt_δs(opt_func, c2, δs)
    n = length(δs)
    xs = Vector{Float64}(undef, n)
    vs = Vector{Float64}(undef, n)
    for (i, δ) in enumerate(δs)
        xs[i], vs[i] = opt_func(c2 + δ, c2)
    end
    return xs, vs
end

const δs = range(0, π, 1001)
const rs_π = opt_δs(opt_angle_π, π, δs)
const rs_2π = opt_δs(opt_angle_2π, 2π, δs)

const prefix = joinpath(@__DIR__, "imgs", "three_pulses")

figure()
plot(rs_π[2], π .+ δs, "C0-", label="\$\\pi\$")
plot(δs, π .+ δs, "C0--")
plot(rs_2π[2], 2π .+ δs, "C1-", label="\$2\\pi\$")
plot(δs, 2π .+ δs, "C1--")
grid()
xlabel("Desired rotation angle (rad)")
ylabel("Total rotation angle (rad)")
legend(ncol=2)
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsData.Fitting: fit_data

function rotation_matrix(α, β, γ)
    cosα = cos(α)
    sinα = sin(α)
    cosβ = cos(β)
    sinβ = sin(β)
    cosγ = cos(γ)
    sinγ = sin(γ)

    sinαsinβ = sinα * sinβ
    cosαsinβ = cosα * sinβ

    return ((cosβ * cosγ, sinαsinβ * cosγ - cosα * sinγ, cosαsinβ * cosγ + sinα * sinγ),
            (cosβ * sinγ, sinαsinβ * sinγ + cosα * cosγ, cosαsinβ * sinγ - sinα * cosγ),
            (-sinβ, sinα * cosβ, cosα * cosβ))
end

function rotate(α, β, γ, x, y, z)
    R = rotation_matrix(α, β, γ)

    R00, R01, R02 = R[1]
    R10, R11, R12 = R[2]
    R20, R21, R22 = R[3]
    return (R00 * x + R01 * y + R02 * z,
            R10 * x + R11 * y + R12 * z,
            R20 * x + R21 * y + R22 * z)
end

translate(x0, y0, z0, x, y, z) = (x0 + x, y0 + y, z0 + z)

# Projecting on a surface at `z0 = -h` and center pixel is xc, yc.
function project(h, xc, yc, x, y, z)
    return (x, y) .* (h / -z) .+ (xc, yc)
end

# params = (α, β, γ, x0, y0, z0)
function image_projection(params, h, xc, yc, p)
    p1 = rotate(params[1], params[2], params[3], p...)
    p2 = translate(params[4], params[5], params[6], p1...)
    return project(h, xc, yc, p2...)
end

# a x + b y = c
# d x + e y = f
function solve_linear(a, b, c, d, e, f)
    denom = (a * e - b * d)
    return (c * e - b * f) / denom, (c * d - a * f) / -denom
end

# (R00 * xb + R01 * yb + R02 * zb + x0) * h + xp * (R20 * xb + R21 * yb + R22 * zb + z0) = 0
# (R10 * xb + R11 * yb + R12 * zb + y0) * h + yp * (R20 * xb + R21 * yb + R22 * zb + z0) = 0

# (R00 * h + xp * R20) * xb + (R01 * h + xp * R21) * yb = -(R02 * h * zb + x0 * h + xp * R22 * zb + xp * z0)
# (R10 * h + yp * R20) * xb + (R11 * h + yp * R21) * yb = -(R12 * h * zb + y0 * h + yp * R22 * zb + yp * z0)

function reverse_projection(params, h, xc, yc, zb, xp, yp)
    R = rotation_matrix(params[1], params[2], params[3])
    x0, y0, z0 = params[4], params[5], params[6]
    xp -= xc
    yp -= yc

    R00, R01, R02 = R[1]
    R10, R11, R12 = R[2]
    R20, R21, R22 = R[3]

    a = R00 * h + xp * R20
    b = R01 * h + xp * R21
    c = -(R02 * h * zb + x0 * h + xp * R22 * zb + xp * z0)
    d = R10 * h + yp * R20
    e = R11 * h + yp * R21
    f = -(R12 * h * zb + y0 * h + yp * R22 * zb + yp * z0)

    return solve_linear(a, b, c, d, e, f)
end

const pixel_coordinates = (
    O = (2085, 1608),
    A = (2094, 1084),
    B = (2283, 1088),
    C = (2294, 562),
    D = (1431, 1037),
    E = (1642, 1043),
    F = (1429, 1158),
    G = (1640, 1163),
    H = (2029, 750),
    I = (2027, 870),
    J = (1411, 858),
    K = (1291, 855),
    L = (959, 888),
    M = (574, 877),
    N = (1575, 1626),
    P = (1573, 1748),
    Q = (978, 1699),
    R = (302, 1689),
    S = (973, 1661),
    T = (305, 1650),
)

const zboard = (
    O = -0.230,
    A = -0.230,
    B = -0.230,
    C = -0.230,
    D = 0.0,
    E = 0.0,
    F = 0.0,
    G = 0.0,
    H = 0.0,
    I = 0.0,
    J = 0.0,
    K = 0.0,
    L = 0.0,
    M = 0.0,
    N = 0.0,
    P = 0.0,
    Q = 0.0,
    R = 0.0,
    S = 0.0,
    T = 0.0,
)

const img_xc = 2016.0
const img_yc = 1512.0
const img_h = 8503.55
@inline function board_pos_for(param, name)
    return reverse_projection(param, img_h, img_xc, img_yc, getfield(zboard, name),
                              getfield(pixel_coordinates, name)...)
end

function dist(p1, p2)
    d = p1 .- p2
    return hypot(d...)
end

function single_model(idx, param)
    if idx == 1
        # O = (0, 0, -0.230)
        return board_pos_for(param, :O)[1] - 0
    elseif idx == -1
        # O = (0, 0, -0.230)
        return board_pos_for(param, :O)[2] - 0
    elseif idx == 2
        # A = (0, -1.850, -0.230)
        return board_pos_for(param, :A)[1] - 0
    elseif idx == -2
        # A = (0, -1.850, -0.230)
        return board_pos_for(param, :A)[2] - -1.850
    elseif idx == 3
        # B = (0.678, -1.850, -0.230)
        return board_pos_for(param, :B)[1] - 0.678
    elseif idx == -3
        # B = (0.678, -1.850, -0.230)
        return board_pos_for(param, :B)[2] - -1.850
    elseif idx == 4
        # C = (0.678, -3.700, -0.230)
        return board_pos_for(param, :C)[1] - 0.678
    elseif idx == -4
        # C = (0.678, -3.700, -0.230)
        return board_pos_for(param, :C)[2] - -3.700
    elseif idx == 5
        # D - E = (-0.738, 0)
        return board_pos_for(param, :D)[1] - board_pos_for(param, :E)[1] - -0.738
    elseif idx == -5
        # D - E = (-0.738, 0)
        return board_pos_for(param, :D)[2] - board_pos_for(param, :E)[2] - 0
    elseif idx == 6
        # D - F = (0, -0.422)
        return board_pos_for(param, :D)[1] - board_pos_for(param, :F)[1] - 0
    elseif idx == -6
        # D - F = (0, -0.422)
        return board_pos_for(param, :D)[2] - board_pos_for(param, :F)[2] - -0.422
    elseif idx == 7
        # F - G = (-0.738, 0)
        return board_pos_for(param, :F)[1] - board_pos_for(param, :G)[1] - -0.738
    elseif idx == -7
        # F - G = (-0.738, 0)
        return board_pos_for(param, :F)[2] - board_pos_for(param, :G)[2] - 0
    elseif idx == 8
        # H - I = (0, -0.422)
        return board_pos_for(param, :H)[1] - board_pos_for(param, :I)[1] - 0
    elseif idx == -8
        # H - I = (0, -0.422)
        return board_pos_for(param, :H)[2] - board_pos_for(param, :I)[2] - -0.422
    elseif idx == 9
        # K - J = (-0.422, 0)
        return board_pos_for(param, :K)[1] - board_pos_for(param, :J)[1] - -0.422
    elseif idx == -9
        # K - J = (-0.422, 0)
        return board_pos_for(param, :K)[2] - board_pos_for(param, :J)[2] - 0
    elseif idx == 10
        # Y(L) - (Y(T) + Y(S) + Y(R) + Y(Q)) / 4 = -2.745
        YL = board_pos_for(param, :L)[2]
        YT = board_pos_for(param, :T)[2]
        YS = board_pos_for(param, :S)[2]
        YR = board_pos_for(param, :R)[2]
        YQ = board_pos_for(param, :Q)[2]
        return YL - (YT + YS + YR + YQ) / 4 - -2.745
    elseif idx == 11
        # Y(M) - Y(L) = 0
        YM = board_pos_for(param, :M)[2]
        YL = board_pos_for(param, :L)[2]
        return YM - YL
    elseif idx == 12
        # Y(T) - Y(S) = 0
        YT = board_pos_for(param, :T)[2]
        YS = board_pos_for(param, :S)[2]
        return YT - YS
    elseif idx == 13
        # Y(R) - Y(Q) = 0
        YR = board_pos_for(param, :R)[2]
        YQ = board_pos_for(param, :Q)[2]
        return YR - YQ
    elseif idx == 14
        # Y(T) + Y(S) - Y(R) - Y(Q) = -0.260
        YT = board_pos_for(param, :T)[2]
        YS = board_pos_for(param, :S)[2]
        YR = board_pos_for(param, :R)[2]
        YQ = board_pos_for(param, :Q)[2]
        return YT + YS - YR - YQ - -0.260
    end
    return 0.0
end

function model(idxs, param)
    return single_model.(idxs, Ref(param))
end

const xindexes = [-9:-1; 1:14]
const yzeros = zeros(length(xindexes))

const fit = fit_data(model, xindexes, yzeros, [0.0, 0.0, 0.0, 0.0, 0.0, -10.0], plotx=false)
function gen_reverse_mapper(param)
    return (x, y) -> reverse_projection(param, img_h, img_xc, img_yc, 0.0, x, y)
end
const reverse_mapper = gen_reverse_mapper(fit.param)

println("Bottom left: $(reverse_mapper(1070, 1549))")
println("Bottom right: $(reverse_mapper(1575, 1626))")
println("Top left: $(reverse_mapper(959, 888))")
println("Top right: $(reverse_mapper(1640, 1163))")
println("Beam position: $(reverse_mapper(1479, 2738))")

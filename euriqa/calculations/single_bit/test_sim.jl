#!/usr/bin/julia

function xy_to_polar(x, y)
    return hypot(x, y), atan(y, x)
end

function compose_sigmas(I, X, Y, Z)
    return [I + Z X - im * Y
            X + im * Y I - Z]
end

function decompose_sigmas(M)
    I = (M[1, 1] + M[2, 2]) / 2
    Z = (M[1, 1] - M[2, 2]) / 2
    X = (M[1, 2] + M[2, 1]) / 2
    Y = (M[1, 2] - M[2, 1]) * im / 2
    return Complex(I), Complex(X), Complex(Y), Complex(Z)
end

function xy_rotation(ψ, θ)
    #   e^(i * (σ_x cosθ + σ_y sinθ) ψ / 2)
    # = cos(ψ / 2) + i * sin(ψ / 2) * (σ_x cosθ + σ_y sinθ)
    # = [cos(ψ / 2)                    i * sin(ψ / 2) (cosθ - i sinθ)
    #    i * sin(ψ / 2) (cosθ + i sinθ)  cos(ψ / 2)]
    sinθ, cosθ = sincos(θ)
    sinψ_2, cosψ_2 = sincos(ψ / 2)
    return [cosψ_2 sinψ_2 * (im * cosθ + sinθ)
            sinψ_2 * (im * cosθ - sinθ) cosψ_2]
end

function z_rotation(ψ)
    #   e^(i * σ_z ψ / 2)
    # = cos(ψ / 2) + i * sin(ψ / 2) * σ_z
    sinψ_2, cosψ_2 = sincos(ψ / 2)
    return compose_sigmas(cosψ_2, 0, 0, im * sinψ_2)
end

function decompose_xy_z(M)
    a0, x0, y0, z0 = decompose_sigmas(M)

    a1 = sqrt(a0^2 - z0^2)
    one_v = one(typeof(a1))
    zero_v = zero(typeof(a1))

    if a1 == 0
        # This has to be the case if the decomposition is possible.
        return ((a0, x0, y0, z0),
                (one_v, zero_v, zero_v, zero_v))
    end
    cosθ = a0 / a1
    isinθ = z0 / a1
    x1 = (a0 * x0 - im * z0 * y0) / a1
    y1 = (a0 * y0 + im * z0 * x0) / a1

    return (a1, x1, y1, zero_v), (cosθ, zero_v, zero_v, isinθ)
end

function equivalent_xy_angle(M)
    σ_xy, = decompose_xy_z(M)
    tanψ_2 = sqrt(abs(σ_xy[2]^2 + σ_xy[3]^2)) / abs(σ_xy[1])
    return atan(tanψ_2) * 2
end

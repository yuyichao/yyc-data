#!/usr/bin/julia

using StaticArrays

function xy_to_polar(x, y)
    return hypot(x, y), atan(y, x)
end

function polar_to_xy(r, θ)
    sinθ, cosθ = sincos(θ)
    return r * cosθ, r * sinθ
end

function compose_sigmas(I, X, Y, Z)
    return SA[I + Z X - im * Y
              X + im * Y I - Z]
end

function decompose_sigmas(M)
    I = (M[1, 1] + M[2, 2]) / 2
    Z = (M[1, 1] - M[2, 2]) / 2
    X = (M[1, 2] + M[2, 1]) / 2
    Y = (M[1, 2] - M[2, 1]) * im / 2
    return Complex(I), Complex(X), Complex(Y), Complex(Z)
end

function compose_xy_rot(ψ, θ)
    #   e^(i * (σ_x cosθ + σ_y sinθ) ψ / 2)
    # = cos(ψ / 2) + i * sin(ψ / 2) * (σ_x cosθ + σ_y sinθ)
    # = [cos(ψ / 2)                    i * sin(ψ / 2) (cosθ - i sinθ)
    #    i * sin(ψ / 2) (cosθ + i sinθ)  cos(ψ / 2)]
    sinθ, cosθ = sincos(θ)
    sinψ_2, cosψ_2 = sincos(ψ / 2)
    return SA[cosψ_2 sinψ_2 * (im * cosθ + sinθ)
              sinψ_2 * (im * cosθ - sinθ) cosψ_2]
end

function decompose_xy_rot((a0, x0, y0, z0)::NTuple{4})
    abs_a0 = abs(a0)
    g = abs_a0 / a0
    x0 = x0 * g
    y0 = y0 * g
    z0 = z0 * g
    a0 = abs_a0
    ψ_2 = atan(sqrt(abs2(x0) + abs2(y0)), a0)
    θ = atan(imag(y0), imag(x0))
    return ψ_2 * 2, θ
end
@inline decompose_xy_rot(M::AbstractMatrix) = decompose_xy_rot(decompose_sigmas(M))

function compose_z_rot(ψ)
    #   e^(i * σ_z ψ / 2)
    # = cos(ψ / 2) + i * sin(ψ / 2) * σ_z
    sinψ_2, cosψ_2 = sincos(ψ / 2)
    return compose_sigmas(cosψ_2, 0, 0, im * sinψ_2)
end

function decompose_z_rot((a0, x0, y0, z0)::NTuple{4})
    return atan(imag(z0 / a0)) * 2
end
@inline decompose_z_rot(M::AbstractMatrix) = decompose_z_rot(decompose_sigmas(M))

function decompose_xy_z((a0, x0, y0, z0)::NTuple{4})
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
@inline decompose_xy_z(M::AbstractMatrix) = decompose_xy_z(decompose_sigmas(M))

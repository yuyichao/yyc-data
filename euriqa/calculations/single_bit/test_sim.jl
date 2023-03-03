#!/usr/bin/julia

function xy_to_polar(x, y)
    return hypot(x, y), atan(y, x)
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

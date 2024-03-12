#!/usr/bin/julia

function get_decay_H(N, dE, Γ)
    dΩ = sqrt(Γ * dE / (2π))
    H = zeros(N + 1, N + 1)
    for i in 1:N
        H[1, i + 1] = dΩ
        H[i + 1, 1] = dΩ
        H[i + 1, i + 1] = (i - (N + 1) / 2) * dE
    end
    return H
end

function draw_x(ψph, k0, dk, xs)
    N = length(ψph)
    return [begin
                v = zero(ComplexF64)
                for (i, ψ) in enumerate(ψph)
                    k = -(k0 + (i - (N + 1) / 2) * dk)
                    v += cis(k * x) * ψ
                end
                v * sqrt(dk)
            end for x in xs]
end

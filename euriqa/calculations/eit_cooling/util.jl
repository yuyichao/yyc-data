#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Trap

using LinearAlgebra

function fill_lambda_H!(H::AbstractMatrix, Γ, Ω₁, Ω₂, Δ₁, Δ₂,
                        ωs::NTuple{N}, ηs::NTuple{N}, nmotions::NTuple{N}) where N
    nmotion = prod(nmotions)
    if size(H) != (nmotion * 3, nmotion * 3)
        error("Hamiltonian dimension is wrong")
    end

    motion_idxs = CartesianIndices(nmotions)
    motion_lidxs = LinearIndices(nmotions)

    @inbounds for idx1 in motion_idxs
        lidx1 = motion_lidxs[idx1]
        n1 = idx1.I .- 1
        for i2 in motion_idxs
            lidx2 = motion_lidxs[idx2]
            n2 = idx2.I .- 1
            # Diagnal
            if n1 == n2
                E_motion = sum(ωs .* n1)
                H[lidx1, lidx1] = Δ₁ + E_motion
                H[lidx1 + nmotion, lidx1 + nmotion] = im * Γ + E_motion
                H[lidx1 + nmotion * 2, lidx1 + nmotion * 2] = Δ₂ + E_motion
            else
                H[lidx1, lidx2] = 0
                H[lidx1 + nmotion, lidx2 + nmotion] = 0
                H[lidx1 + nmotion * 2, lidx2 + nmotion * 2] = 0
            end
            # 1-3, uncoupled
            H[lidx1, lidx2 + nmotion * 2] = 0
            H[lidx1 + nmotion * 2, lidx2] = 0

            # 1-2 and 2-3
            M = prod(Trap.sideband_with_phase.(n1, n2, ηs))
            H[lidx1, lidx2 + nmotion] .= Ω₁ * M
            H[lidx1 + nmotion, lidx2 + nmotion * 2] .= Ω₂ * M
        end
    end
end

# function Λ_scatter(Γ, Ω₁, Ω₂, Δ₁, Δ₂)
#     M = [    Δ₁   Ω₁ / 2        0
#          Ω₁ / 2   im * Γ   Ω₂ / 2
#               0   Ω₂ / 2       Δ₂]
#     E = eigen(M)
#     vecs = E.vectors
#     vals = E.values

#     idx = 0
#     max_prob = 0.0
#     for i in 1:3
#         prob = abs2(vecs[1, i])
#         if prob >= max_prob
#             max_prob = prob
#             idx = i
#         end
#     end
#     return imag(vals[idx])
# end

# function Λ_scatter(Γ, Ω₁, Ω₂, Δ₁, Δ₂)
#     M = [    Δ₁   Ω₁ / 2        0
#          Ω₁ / 2   im * Γ   Ω₂ / 2
#               0   Ω₂ / 2       Δ₂]
#     E = eigen(M)
#     vecs = E.vectors
#     vals = E.values

#     idx = 0
#     max_prob = 0.0
#     for i in 1:3
#         prob = abs2(vecs[2, i])
#         if prob >= max_prob
#             max_prob = prob
#             idx = i
#         end
#     end
#     popat!(vals, idx)
#     g1 = imag(vals[1])
#     g2 = imag(vals[2])
#     return g1 * g2 / (g1 + g2)
# end

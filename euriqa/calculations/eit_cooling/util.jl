#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Trap

using LinearAlgebra

mutable struct SidebandCache{T}
    η::T
    const values::Vector{Complex{T}}
    function SidebandCache{T}() where T
        return new(0, Complex{T}[])
    end
end

# Have these inlined somehow cause the parent function to bloat
@noinline function _clear_sideband!(cache::SidebandCache, η)
    empty!(cache.values)
    cache.η = η
end
@noinline function _resize_sideband!(values::Vector, new_nvalues, nvalues)
    resize!(values, new_nvalues)
    @inbounds values[nvalues + 1:new_nvalues] .= NaN
    return
end
@noinline function _fill_sideband_value!(values::Vector{CT}, idx, n1, n2, η) where CT
    val = CT(Trap.sideband_with_phase(n1, n2, η))
    @inbounds values[idx] = val
    return val
end

@inline function get_sideband(cache::SidebandCache, n1, n2, η)
    if n1 < n2
        n1, n2 = n2, n1
    end
    n1_offset = (n1 * (n1 + 1)) >> 1
    idx = n1_offset + 1 + n2

    if cache.η != η
        _clear_sideband!(cache, η)
    end

    values = cache.values
    nvalues = length(values)
    if idx > nvalues
        _resize_sideband!(values, n1_offset + n1 + 1, nvalues)
    end
    cached_val = @inbounds values[idx]
    if isnan(real(cached_val))
        cached_val = _fill_sideband_value!(values, idx, n1, n2, η)
    end
    return cached_val
end

@noinline function fill_upto!(cache::SidebandCache, n, η)
    if cache.η != η
        if η == 0
            return
        end
        _clear_sideband!(cache, η)
    end
    values = cache.values
    nvalues = length(values)
    new_nvalues = (n + 1) * (n + 2) >> 1 + 1
    if new_nvalues > nvalues
        _resize_sideband!(values, new_nvalues, nvalues)
    end
    idx = 0
    @inbounds for n1 in 0:n
        for n2 in 0:n1
            idx += 1
            cached_val = values[idx]
            if !isnan(real(cached_val))
                continue
            end
            values[idx] = Trap.sideband_with_phase(n1, n2, η)
        end
    end
end

@inline function get_sideband_nocheck(cache::SidebandCache{T}, n1, n2, η) where T
    if n1 < n2
        n1, n2 = n2, n1
    end
    if η == 0
        return Complex{T}(n1 == n2 ? 1 : 0)
    end
    return @inbounds cache.values[(n1 * (n1 + 1)) >> 1 + 1 + n2]
end

struct Builder{T,N}
    sideband_caches::NTuple{N,SidebandCache{T}}
    function Builder{T,N}() where {T,N}
        return new(ntuple(x->SidebandCache{T}(), Val(N)))
    end
end

function fill_lambda_H!(builder::Builder{T,N}, H::AbstractMatrix, Γ, Ω₁, Ω₂, Δ₁, Δ₂,
                        ωs::NTuple{N}, ηs::NTuple{N}, nmotions::NTuple{N}) where {T,N}
    nmotion = prod(nmotions)
    if size(H) != (nmotion * 3, nmotion * 3)
        error("Hamiltonian dimension is wrong")
    end

    motion_idxs = CartesianIndices(nmotions)
    motion_lidxs = LinearIndices(nmotions)

    sideband_caches = builder.sideband_caches

    fill_upto!.(sideband_caches, nmotions .- 1, ηs)

    @inbounds for idx2 in motion_idxs
        lidx2 = motion_lidxs[idx2]
        n2 = idx2.I .- 1
        for idx1 in motion_idxs
            lidx1 = motion_lidxs[idx1]
            n1 = idx1.I .- 1
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
            M = prod(get_sideband_nocheck.(sideband_caches, n1, n2, ηs))
            H[lidx1, lidx2 + nmotion] = Ω₁ * M
            H[lidx1 + nmotion, lidx2 + nmotion * 2] = Ω₂ * M
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

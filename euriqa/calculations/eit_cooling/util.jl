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
                H[lidx1 + nmotion, lidx1 + nmotion] = im * Γ / 2 + E_motion
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

function fill_scatter_branch!(builder::Builder{T,N}, B::AbstractMatrix, ηs::NTuple{N},
                              nmotions::NTuple{N}) where {T,N}
    nmotion = prod(nmotions)
    if size(B) != (nmotion, nmotion)
        error("Branching matrix dimension is wrong")
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

            B[lidx1, lidx2] = prod(get_sideband_nocheck.(sideband_caches, n1, n2, ηs))
        end
    end
end

@kwdef struct DecayParams{T,N}
    p::T
    amp::NTuple{2,T} # Amplitude of 1/2 ground state
    ηs::NTuple{N,T}
end

@kwdef struct SystemParams{T,N}
    Γ::T
    Ω₁::T
    Ω₂::T
    Δ₁::T
    Δ₂::T
    ωs::NTuple{N,T}
    ηs::NTuple{N,T}
    nmotions::NTuple{N,Int}
    decay_branch::Vector{DecayParams{T,N}} = DecayParams{T,N}[]
end

struct RateMatrices{T,N}
    R::Matrix{T}
    U::Matrix{T}
end

function build_rate_matrices(builder::Builder{T,N},
                             params::SystemParams{T,N}) where {T,N}
    nmotion = prod(params.nmotions)
    H = Matrix{Complex{T}}(undef, nmotion * 3, nmotion * 3)
    fill_lambda_H!(builder, H, params.Γ, params.Ω₁, params.Ω₂, params.Δ₁, params.Δ₂,
                   params.ωs, params.ηs, params.nmotions)
    F = eigen!(H)

    U = F.vectors # Unitary from new basis to old basis
    Ud = U' # Unitary from old basis to new basis

    rates = imag.(F.values) .* 2

    # Rates
    R = zeros(T, nmotion * 3, nmotion * 3)
    B = Matrix{Complex{T}}(undef, nmotion, nmotion)
    tmp = Matrix{Complex{T}}(undef, nmotion * 3, nmotion)

    for decay in params.decay_branch
        fill_scatter_branch!(builder, B, decay.ηs, params.nmotions)
        mul!(tmp, @view(Ud[:, nmotion + 1:nmotion * 2]), B)
        tmp .*= rates
        if decay.amp[1] > 0
            mul!(H, tmp, @view(U[1:nmotion, :]), decay.amp[1], 0)
            if decay.amp[2] > 0
                mul!(H, tmp, @view(U[nmotion * 2 + 1:nmotion * 3, :]),
                     decay.amp[2], 1)
            end
        else
            mul!(H, tmp, @view(U[nmotion * 2 + 1:nmotion * 3, :]),
                 decay.amp[2], 0)
        end
        R .+= decay.p .* abs2.(H)
    end

    for i in 1:nmotion * 3
        R[i, i] -= rates[i]
    end

    return RateMatrices{T,N}(R, abs2.(U))
end

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

@kwdef struct DecayParams{T,N}
    p::T
    amp::NTuple{2,T} # Amplitude of 1/2 ground state
    ηs::NTuple{N,T}
end

@kwdef mutable struct SystemParams{T,N}
    Γ::T
    Ω₁::T
    Ω₂::T
    Δ₁::T
    Δ₂::T
    ωs::NTuple{N,T}
    ηs::NTuple{N,T}
    const nmotions::NTuple{N,Int}
    const decay_branch::Vector{DecayParams{T,N}} = DecayParams{T,N}[]
end

struct Builder{T,N}
    params::SystemParams{T,N}
    H::Matrix{Complex{T}}
    B::Matrix{Complex{T}}
    tmp::Matrix{Complex{T}}
    tmp2::Vector{T}
    sideband_caches::NTuple{N,SidebandCache{T}}
    function Builder(params::SystemParams{T,N}) where {T,N}
        nmotion = prod(params.nmotions, init=1)
        return new{T,N}(params,
                        Matrix{Complex{T}}(undef, nmotion * 3, nmotion * 3),
                        Matrix{Complex{T}}(undef, nmotion, nmotion),
                        Matrix{Complex{T}}(undef, nmotion, nmotion * 3),
                        Vector{T}(undef, nmotion * 3),
                        ntuple(x->SidebandCache{T}(), Val(N)))
    end
end

function _fill_H!(builder::Builder{T,N}) where {T,N}
    H = builder.H
    params = builder.params
    Γ = params.Γ
    Ω₁ = params.Ω₁
    Ω₂ = params.Ω₂
    Δ₁ = params.Δ₁
    Δ₂ = params.Δ₂
    ωs = params.ωs
    ηs = params.ηs
    nmotions = params.nmotions
    nmotion = prod(nmotions, init=1)

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
                E_motion = sum(ωs .* n1, init=zero(T))
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
            M = prod(get_sideband_nocheck.(sideband_caches, n1, n2, ηs),
                     init=one(Complex{T}))
            H[lidx1, lidx2 + nmotion] = Ω₁ * M / 2
            H[lidx1 + nmotion, lidx2] = Ω₁ * M / 2
            H[lidx1 + nmotion, lidx2 + nmotion * 2] = Ω₂ * M / 2
            H[lidx1 + nmotion * 2, lidx2 + nmotion] = Ω₂ * M / 2
        end
    end
end

function _fill_scatter_branch!(builder::Builder{T,N}, ηs::NTuple{N}) where {T,N}
    B = builder.B
    nmotions = builder.params.nmotions

    nmotion = prod(nmotions, init=1)
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

            B[lidx1, lidx2] = prod(get_sideband_nocheck.(sideband_caches, n1, n2, ηs),
                                   init=one(Complex{T}))
        end
    end
end

struct RateMatrices{T,N}
    R::Matrix{T}
    U::Matrix{T}
    rates::Vector{T}
end

function propagate(cb::F, p0, rates::RateMatrices, ts) where F
    cur_t = Ref(zero(ts[1]))
    _p = Ref(rates.U \ p0)

    @inline function compute_result(t)
        dt = t - cur_t[]
        if dt != 0
            _p[] = LinearAlgebra.exp!(dt .* rates.R) * _p[]
            cur_t[] = t
        end
        return cb(t, rates.U * _p[], rates)
    end
    return compute_result.(ts)
end

function build_rate_matrices(builder::Builder{T,N}) where {T,N}
    _fill_H!(builder)

    params = builder.params
    H = builder.H
    F = eigen!(H)

    U = F.vectors # From new basis to old basis
    Ud = inv(U) # From old basis to new basis

    rates = imag.(F.values) .* 2

    # Rates
    nmotion = prod(params.nmotions, init=1)
    R = zeros(T, nmotion * 3, nmotion * 3)
    B = builder.B
    tmp = builder.tmp
    tmp2 = builder.tmp2

    for decay in params.decay_branch
        _fill_scatter_branch!(builder, decay.ηs)
        mul!(tmp, B, @view(U[nmotion + 1:nmotion * 2, :]))
        if decay.amp[1] > 0
            mul!(H, @view(Ud[:, 1:nmotion]), tmp, decay.amp[1], false)
            if decay.amp[2] > 0
                mul!(H, @view(Ud[:, nmotion * 2 + 1:nmotion * 3]),
                     tmp, decay.amp[2], true)
            end
        else
            mul!(H, @view(Ud[:, nmotion * 2 + 1:nmotion * 3]),
                 tmp, decay.amp[2], false)
        end

        @inbounds for j in 1:nmotion * 3
            s = zero(T)
            r = rates[j]
            @simd for i in 1:nmotion * 3
                p = abs2(H[i, j])
                tmp2[i] = p
                s += p
            end
            if s == 0
                continue
            end
            scale = r / s * decay.p
            @simd for i in 1:nmotion * 3
                R[i, j] = muladd(scale, tmp2[i], R[i, j])
            end
            R[j, j] -= r
        end
    end

    return RateMatrices{T,N}(R, abs2.(U), rates)
end

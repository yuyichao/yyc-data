#!/usr/bin/julia

module Utils

using SpecialFunctions

# LGPLv3 implementation from libstdc++

function poly_laguerre_large_n(n::Integer, α, x::Tp) where Tp
    a::Tp = -n
    b::Tp = α + 1
    η::Tp = 2b - 4a
    cos²th = x / η
    sin²th = 1 - cos²th
    costh = @fastmath sqrt(cos²th)
    th = @fastmath acos(costh)
    pre_h = (Tp(π / 2)^2) * η * η * cos²th * sin²th
    lg_b = @fastmath lgamma(Tp(n + b))
    lnfact = @fastmath lgamma(Tp(n + 1))
    pre_term1 = @fastmath Tp(0.5) * (1 - b) * log(Tp(0.25) * x * η)
    pre_term2 = @fastmath Tp(0.25) * log(pre_h)
    lnpre = lg_b - lnfact + Tp(0.5) * x + pre_term1 - pre_term2

    th2 = 2 * th
    sin2th = @fastmath 2 * costh * sqrt(sin²th)

    # From libstdc++
    ser_term1 = 0 # @fastmath sinpi(a)
    # This might be off by a minus sign or sth like that.
    # Evaluating at `10000001, 10, 1.2` gives the wrong result
    ser_term2 = @fastmath sin(Tp(0.25) * η * (th2 - sin2th) + Tp(π / 4))

    ser::Tp = ser_term1 + ser_term2
    return @fastmath exp(lnpre) * ser
end

function poly_laguerre_hyperg(n::Integer, α, x::Tp) where Tp
    b::Tp = Tp(α) + 1
    mx = -x
    tc_sgn::Tp = x < 0 ? 1 : ((n % 2 == 1) ? -1 : 1)
    # Get |x|^n/n!
    tc::Tp = 1
    ax = abs(x)
    for k in 1:n
        tc *= ax / k
    end
    term::Tp = tc * tc_sgn
    _sum::Tp = term
    for k in (n - 1):-1:0
        term *= ((b + Tp(k)) / Tp(n - k)) * Tp(k + 1) / mx
        _sum += term
    end
    return _sum
end

function poly_laguerre_recursion(n::Integer, α, x::Tp) where Tp
    # Compute l_0.
    l_0::Tp = 1
    n == 0 && return l_0

    # Compute l_1^alpha.
    l_1::Tp = -x + 1 + α
    n == 1 && return l_1

    # Compute l_n^alpha by recursion on n.
    l_n′::Tp = l_0
    l_n::Tp = l_1
    b::Tp = α - 1
    a::Tp = b - x
    @fastmath for nn in 2:n
        fnn = Tp(nn)
        l1 = muladd(a, l_n, -b * l_n′)
        l2 = muladd(2, l_n, -l_n′)
        l_n, l_n′ = l1 / fnn + l2, l_n
    end
    return l_n
end

function (genlaguerre(n::Integer, α, x::Tp)::Tp) where Tp<:AbstractFloat
    if x < 0
        throw(DomainError())
    elseif isnan(α) || isnan(x)
        # Return NaN on NaN input.
        return NaN
    elseif n == 0
        return 1
    elseif n == 1
        return Tp(1) + Tp(α) - x
    elseif x == 0
        prod::Tp = α + 1
        for k in 2:n
            prod *= Tp(α + k) / Tp(k)
        end
        return prod
    elseif n > 10000000 && α > -1 && x < 2 * (α + 1) + 4n
        return poly_laguerre_large_n(n, α, x)
    elseif α >= 0 || (x > 0 && α < -(n + 1))
        return poly_laguerre_recursion(n, α, x)
    else
        return poly_laguerre_hyperg(n, α, x)
    end
end
genlaguerre(n::Integer, α, x) = genlaguerre(n, α, float(x))

function binomial_estimate(x, n, z::T=1.0) where T<:AbstractFloat
    if n <= 0
        return T(0.5), T(0.5)
    end
    p = T(x / n)
    z² = z^2
    z²n = z² / n
    p′::T = (p + z²n / 2) / (1 + z²n)
    unc::T = sqrt(p * (1 - p) / n + z² / 4 / n^2) / (1 + z²n)
    return p′, unc
end

function binomial_interval(x, n, z::T=1.0) where T<:AbstractFloat
    p, unc = binomial_estimate(x, n, z)
    return p - unc, p + unc
end

linspace(a, b, n) = range(a, stop=b, length=n)

using Random

const _interactive = Ref(true)

thread_rng() = Random.GLOBAL_RNG

_nbits(::Type{Bool}) = 1
_nbits(::Type{T}) where T = sizeof(T) * 8

# Copied from randsubseq! from Base
@inline function __rand_setbits(r::AbstractRNG, ::Type{T}, L) where T
    n = sizeof(T) * 8
    S = zero(T)
    # Skip through A, in order, from each element i to the next element i+s
    # included in S. The probability that the next included element is
    # s==k (k > 0) is (1-p)^(k-1) * p, and hence the probability (CDF) that
    # s is in {1,...,k} is 1-(1-p)^k = F(k).   Thus, we can draw the skip s
    # from this probability distribution via the discrete inverse-transform
    # method: s = ceil(F^{-1}(u)) where u = rand(), which is simply
    # s = ceil(log(rand()) / log1p(-p)).
    # -log(rand()) is an exponential variate, so can use randexp().
    i = 0
    while true
        s = randexp(r) * L
        s >= n - i && return S # compare before ceil to avoid overflow
        i += unsafe_trunc(Int, ceil(s))
        S |= one(T) << (i - 1)
    end
    # [This algorithm is similar in spirit to, but much simpler than,
    #  the one by Vitter for a related problem in "Faster methods for
    #  random sampling," Comm. ACM Magazine 7, 703-718 (1984).]
end

@inline function _rand_setbits(r::AbstractRNG, ::Type{T}, p::Real) where T
    L = -1 / log1p(-p) # L > 0
    return __rand_setbits(r, T, L)
end

mutable struct RandSetBits{T<:Union{Bool,Base.BitInteger},R<:AbstractRNG,P<:Real}
    state::UInt128
    ele_left::Int
    # if p_L >= 0, it is the original probability
    # if p_L < 0, it's 1 / log1p(-p) used for computing the amount to jump ahead.
    const p_L::P
    const rng::R
    function RandSetBits{T}(r::R, p::P) where {T<:Union{Bool,Base.BitInteger},R<:AbstractRNG,P<:Real}
        0 <= p <= 1 || throw(ArgumentError("probability $p not in [0,1]"))
        if 0 < p <= 0.17 # empirical threshold for trivial O(n) algorithm to be better
            p = 1 / log1p(-p) # L > 0
        end
        return new{T,R,P}(0, 0, p, r)
    end
    RandSetBits{T}(p::Real) where T<:Union{Bool,Base.BitInteger} =
        RandSetBits{T}(Random.default_rng(), p)
end

function _new_state(sb::RandSetBits)
    Te = UInt128
    n = sizeof(Te) * 8
    S = zero(Te)
    r = sb.rng
    p = sb.p_L
    @inline if p >= 0
        if p == 0
            return S
        elseif p == 1
            return ~S
        end
        for i = 1:n
            S |= Te(rand(r) <= p) << (i - 1)
        end
    else
        S = __rand_setbits(r, Te, -p)
    end
    return S
end

@inline function Random.rand(sb::RandSetBits{T}) where T
    Te = UInt128
    @inline if sizeof(T) == sizeof(Te)
        return _new_state(sb) % T
    end
    n = sizeof(Te) * 8
    state = sb.state
    if sb.ele_left == 0
        sb.ele_left = n ÷ _nbits(T)
        state = _new_state(sb)
    end
    sb.ele_left -= 1
    res = state % T
    sb.state = state >> _nbits(T)
    return res
end

function rand_setbits(r::AbstractRNG, ::Type{T}, p::Real) where {T <: Union{Bool,Base.BitInteger}}
    0 <= p <= 1 || throw(ArgumentError("probability $p not in [0,1]"))
    @inline if T === Bool
        return rand(r) <= p
    end
    n = sizeof(T) * 8
    p == 1 && return ~zero(T)
    S = zero(T)
    p == 0 && return S
    @inline if p > 0.17 # empirical threshold for trivial O(n) algorithm to be better
        for i = 1:n
            S |= T(rand(r) <= p) << (i - 1)
        end
        return S
    else
        return _rand_setbits(r, T, p)
    end
end
@inline rand_setbits(::Type{T}, p::Real) where T =
    rand_setbits(Random.default_rng(), T, p)


mutable struct RandDepol{T<:Union{Bool,Base.BitInteger},R<:AbstractRNG,P<:Real}
    state_x::UInt128
    state_z::UInt128
    ele_left::Int
    err_sel_left::Int
    err_sel_state::UInt64
    # if p_L >= 0, it is the original probability
    # if p_L < 0, it's 1 / log1p(-p) used for computing the amount to jump ahead.
    const p_L::P
    const rng::R

    function RandDepol{T}(r::R, p::P) where {T<:Union{Bool,Base.BitInteger},R<:AbstractRNG,P<:Real}
        0 <= p <= 1 || throw(ArgumentError("probability $p not in [0,1]"))
        if 0 < p <= 0.17 # empirical threshold for trivial O(n) algorithm to be better
            p = 1 / log1p(-p) # L > 0
        end
        return new{T,R,P}(0, 0, 0, 0, 0, p, r)
    end
    RandDepol{T}(p::Real) where T<:Union{Bool,Base.BitInteger} =
        RandDepol{T}(Random.default_rng(), p)
end

@inline function _rand_err_sel(rd::RandDepol)
    @inline if rd.err_sel_left == 0
        rd.err_sel_state = rand(rd.rng, zero(UInt64):UInt64(UInt64(3)^36 - 1))
        rd.err_sel_left = 36
    end
    rd.err_sel_state, res = divrem(rd.err_sel_state, 3)
    rd.err_sel_left -= 1
    return res % Int
end

function _new_state(rd::RandDepol)
    Te = UInt128
    n = sizeof(Te) * 8
    S1 = zero(Te)
    S2 = zero(Te)
    r = rd.rng
    p = rd.p_L
    @inline if p >= 0
        if p == 0
            return S1, S2
        end
        for i = 1:n
            v = rand(r)
            xerr = v <= p * 2 / 3
            zerr = p / 3 < v <= p
            S1 |= Te(xerr) << (i - 1)
            S2 |= Te(zerr) << (i - 1)
        end
    else
        i = 0
        while true
            s = randexp(r) * -p
            s >= n - i && return S1, S2 # compare before ceil to avoid overflow
            i += unsafe_trunc(Int, ceil(s))
            err_type = _rand_err_sel(rd)
            S1 |= Te(err_type <= 1) << (i - 1)
            S2 |= Te(err_type >= 1) << (i - 1)
        end
    end
    return S1, S2
end

@inline function Random.rand(rd::RandDepol{T}) where T
    Te = UInt128
    @inline if sizeof(T) == sizeof(Te)
        x, z = _new_state(rd)
        return x % T, z % T
    end
    n = sizeof(Te) * 8
    state_x = rd.state_x
    state_z = rd.state_z
    if rd.ele_left == 0
        rd.ele_left = n ÷ _nbits(T)
        state_x, state_z = _new_state(rd)
    end
    rd.ele_left -= 1
    res_x = state_x % T
    res_z = state_z % T
    rd.state_x = state_x >> _nbits(T)
    rd.state_z = state_z >> _nbits(T)
    return res_x, res_z
end


mutable struct Rand2QDepol{T<:Union{Bool,Base.BitInteger},R<:AbstractRNG,P<:Real}
    state_x1::UInt128
    state_z1::UInt128
    state_x2::UInt128
    state_z2::UInt128
    ele_left::Int
    err_sel_left::Int
    err_sel_state::UInt64
    # if p_L >= 0, it is the original probability
    # if p_L < 0, it's 1 / log1p(-p) used for computing the amount to jump ahead.
    const p_L::P
    const rng::R

    function Rand2QDepol{T}(r::R, p::P) where {T<:Union{Bool,Base.BitInteger},R<:AbstractRNG,P<:Real}
        0 <= p <= 1 || throw(ArgumentError("probability $p not in [0,1]"))
        if 0 < p <= 0.17 # empirical threshold for trivial O(n) algorithm to be better
            p = 1 / log1p(-p) # L > 0
        end
        return new{T,R,P}(0, 0, 0, 0, 0, 0, 0, p, r)
    end
    Rand2QDepol{T}(p::Real) where T<:Union{Bool,Base.BitInteger} =
        Rand2QDepol{T}(Random.default_rng(), p)
end

@inline function _rand_err_sel(rd::Rand2QDepol)
    @inline if rd.err_sel_left == 0
        rd.err_sel_state = rand(rd.rng, zero(UInt64):UInt64(UInt64(15)^14 - 1))
        rd.err_sel_left = 14
    end
    rd.err_sel_state, res = divrem(rd.err_sel_state, 15)
    rd.err_sel_left -= 1
    return res % Int
end

function _new_state(rd::Rand2QDepol)
    Te = UInt128
    n = sizeof(Te) * 8
    S1 = zero(Te)
    S2 = zero(Te)
    S3 = zero(Te)
    S4 = zero(Te)
    r = rd.rng
    p = rd.p_L
    @inline if p >= 0
        if p == 0
            return S1, S2, S3, S4
        end
        for i = 1:n
            v = rand(r)
            if v >= p
                continue
            end
            v = unsafe_trunc(Int, floor(v * (15 / p))) + 1
            S1 |= Te(v & 1) << (i - 1)
            S2 |= Te((v >> 1) & 1) << (i - 1)
            S3 |= Te((v >> 2) & 1) << (i - 1)
            S4 |= Te(v >> 3) << (i - 1)
        end
    else
        i = 0
        while true
            s = randexp(r) * -p
            s >= n - i && return S1, S2, S3, S4 # compare before ceil to avoid overflow
            i += unsafe_trunc(Int, ceil(s))
            v = _rand_err_sel(rd) + 1
            S1 |= Te(v & 1) << (i - 1)
            S2 |= Te((v >> 1) & 1) << (i - 1)
            S3 |= Te((v >> 2) & 1) << (i - 1)
            S4 |= Te(v >> 3) << (i - 1)
        end
    end
    return S1, S2, S3, S4
end

@inline function Random.rand(rd::Rand2QDepol{T}) where T
    Te = UInt128
    @inline if sizeof(T) == sizeof(Te)
        x1, z1, x2, z2 = _new_state(rd)
        return x1 % T, z1 % T, x2 % T, z2 % T
    end
    n = sizeof(Te) * 8
    state_x1 = rd.state_x1
    state_z1 = rd.state_z1
    state_x2 = rd.state_x2
    state_z2 = rd.state_z2
    if rd.ele_left == 0
        rd.ele_left = n ÷ _nbits(T)
        state_x1, state_z1, state_x2, state_z2 = _new_state(rd)
    end
    rd.ele_left -= 1
    res_x1 = state_x1 % T
    res_z1 = state_z1 % T
    res_x2 = state_x2 % T
    res_z2 = state_z2 % T
    rd.state_x1 = state_x1 >> _nbits(T)
    rd.state_z1 = state_z1 >> _nbits(T)
    rd.state_x2 = state_x2 >> _nbits(T)
    rd.state_z2 = state_z2 >> _nbits(T)
    return res_x1, res_z1, res_x2, res_z2
end

function __init__()
    interactive_str = get(ENV, "NACS_INTERACT", "true")
    if interactive_str == "false" || interactive_str == "0"
        _interactive[] = false
    end
end

@inline interactive() = _interactive[]

end

#!/usr/bin/julia

module Clifford

using Random
using ..Utils

##
# Implement simulation of clifford circuit using the stabilizer state description

##
# Gates
#
# Gate operations are implemented via their effect on pauli operators.
# The pauli operator is described using two bits, `x` and `z`:
#     x=0, z=0: I
#     x=1, z=0: X
#     x=1, z=1: Y
#     x=0, z=1: Z
# and the `r` carries the sign of the operator (e.g. `x=0`, `z=1`, `r=1` represents
# the `-Z` operator).
#
# For efficient implementation of propagating multiple pauli strings/stabilizers
# at the same time, we allow each of the bits to be packed into an integer
# and we use bitwise operation to implement the gate in parallel.

"""
    nbits(::Type)

Return number of bits packed in the element type.
"""
nbits(::Type{T}) where T = Utils._nbits(eltype(T))

"""
    count_ones(v)

Like `Base.count_ones` but works for `Bool` as well.
"""
count_ones(v::Bool) = Int(v)
count_ones(v) = Base.count_ones(v)

## Single qubit gates
# For these gates, we also support a variety of operations on them for convenience
# including multiplication and inverse.
abstract type Clifford1Q end

struct Composite1Q{G1<:Clifford1Q,G2<:Clifford1Q} <: Clifford1Q
    g1::G1
    g2::G2
end
Base.@propagate_inbounds @inline function apply!(gate::Composite1Q, xas, zas, rs)
    apply!(gate.g1, xas, zas, rs)
    apply!(gate.g2, xas, zas, rs)
    return
end
function Base.:*(g1::Clifford1Q, g2::Clifford1Q)
    return Composite1Q(g1, g2)
end
@inline Base.inv(g::Composite1Q) = Composite1Q(inv(g.g2), inv(g.g1))

# Identity
struct IGate <: Clifford1Q
end
@inline function apply!(::IGate, xas, zas, rs)
    return
end
@inline Base.inv(::IGate) = IGate()
Base.:*(g1::Clifford1Q, ::IGate) = g1
Base.:*(::IGate, g2::Clifford1Q) = g2
Base.:*(::IGate, ::IGate) = IGate()

# Hadamard
struct HGate <: Clifford1Q
end
# I -> I, X -> Z, Y -> -Y, Z -> X
Base.@propagate_inbounds @inline function apply!(::HGate, xas, zas, rs)
    xa = xas[]
    za = zas[]
    r = rs[]

    rs[] = r ⊻ (xa & za)
    xas[] = za
    zas[] = xa
    return
end
@inline Base.inv(::HGate) = HGate()
Base.:*(::HGate, ::HGate) = IGate()

# X
struct XGate <: Clifford1Q
end
# I -> I, X -> X, Y -> -Y, Z -> -Z
Base.@propagate_inbounds @inline function apply!(::XGate, xas, zas, rs)
    za = zas[]
    r = rs[]

    rs[] = r ⊻ za
    return
end
@inline Base.inv(::XGate) = XGate()
Base.:*(::XGate, ::XGate) = IGate()

# Y
struct YGate <: Clifford1Q
end
# I -> I, X -> -X, Y -> Y, Z -> -Z
Base.@propagate_inbounds @inline function apply!(::YGate, xas, zas, rs)
    xa = xas[]
    za = zas[]
    r = rs[]

    rs[] = r ⊻ (xa ⊻ za)
    return
end
@inline Base.inv(::YGate) = YGate()
Base.:*(::YGate, ::YGate) = IGate()

# Z
struct ZGate <: Clifford1Q
end
# I -> I, X -> -X, Y -> -Y, Z -> Z
Base.@propagate_inbounds @inline function apply!(::ZGate, xas, zas, rs)
    xa = xas[]
    r = rs[]

    rs[] = r ⊻ xa
    return
end
@inline Base.inv(::ZGate) = ZGate()
Base.:*(::ZGate, ::ZGate) = IGate()

# Phase gate (sqrt(Z))
struct SGate <: Clifford1Q
end
# I -> I, X -> Y, Y -> -X, Z -> Z
Base.@propagate_inbounds @inline function apply!(::SGate, xas, zas, rs)
    xa = xas[]
    za = zas[]
    r = rs[]

    rs[] = r ⊻ (xa & za)
    zas[] = za ⊻ xa
    return
end
# Inverse of phase gate
struct ISGate <: Clifford1Q
end
# I -> I, X -> -Y, Y -> X, Z -> Z
Base.@propagate_inbounds @inline function apply!(::ISGate, xas, zas, rs)
    xa = xas[]
    za = zas[]
    r = rs[]

    rs[] = r ⊻ (xa & ~za)
    zas[] = za ⊻ xa
    return
end
@inline Base.inv(::SGate) = ISGate()
@inline Base.inv(::ISGate) = SGate()
Base.:*(::SGate, ::SGate) = ZGate()
Base.:*(::ISGate, ::ISGate) = ZGate()

Base.:*(::SGate, ::ZGate) = ISGate()
Base.:*(::ZGate, ::SGate) = ISGate()
Base.:*(::ISGate, ::ZGate) = SGate()
Base.:*(::ZGate, ::ISGate) = SGate()

# sqrt(X)
struct SXGate <: Clifford1Q
end
# I -> I, X -> X, Y -> Z, Z -> -Y
Base.@propagate_inbounds @inline function apply!(::SXGate, xas, zas, rs)
    xa = xas[]
    za = zas[]
    r = rs[]

    rs[] = r ⊻ (~xa & za)
    xas[] = xa ⊻ za
    return
end
# Inverse of sqrt(X)
struct ISXGate <: Clifford1Q
end
# I -> I, X -> X, Y -> -Z, Z -> Y
Base.@propagate_inbounds @inline function apply!(::ISXGate, xas, zas, rs)
    xa = xas[]
    za = zas[]
    r = rs[]

    rs[] = r ⊻ (xa & za)
    xas[] = xa ⊻ za
    return
end
@inline Base.inv(::SXGate) = ISXGate()
@inline Base.inv(::ISXGate) = SXGate()
Base.:*(::SXGate, ::SXGate) = XGate()
Base.:*(::ISXGate, ::ISXGate) = XGate()

Base.:*(::SXGate, ::XGate) = ISXGate()
Base.:*(::XGate, ::SXGate) = ISXGate()
Base.:*(::ISXGate, ::XGate) = SXGate()
Base.:*(::XGate, ::ISXGate) = SXGate()

## Two qubit gates
abstract type Clifford2Q end

struct CNOTGate <: Clifford2Q
end
Base.@propagate_inbounds @inline function apply!(::CNOTGate, xas, zas, xbs, zbs, rs)
    xa = xas[]
    xb = xbs[]
    za = zas[]
    zb = zbs[]
    r = rs[]

    rs[] = r ⊻ (xa & zb & ~(xb ⊻ za))
    xbs[] = xb ⊻ xa
    zas[] = za ⊻ zb
    return
end

# Gate operation function that takes and returns the values
# rather than mutating them in place.
function apply(gate::Clifford1Q, xa, za, r)
    xas = Ref(xa)
    zas = Ref(za)
    rs = Ref(r)
    apply!(gate, xas, zas, rs)
    return xas[], zas[], rs[]
end

function apply(gate::Clifford2Q, xa, za, xb, zb, r)
    xas = Ref(xa)
    zas = Ref(za)
    xbs = Ref(xb)
    zbs = Ref(zb)
    rs = Ref(r)
    apply!(gate, xas, zas, xbs, zbs, rs)
    return xas[], zas[], xbs[], zbs[], rs[]
end

@noinline throw_qubit_bounderror(a) =
    throw(ArgumentError("Qubit index $a out of bound."))
@inline function check_qubit_bound(n, a)
    if ((a - 1) % UInt) >= (n % UInt)
        throw_qubit_bounderror(a)
    end
end

# A single or a set of independent pauli strings
struct PauliString{T<:Union{Bool,Base.BitInteger}}
    n::Int
    xs::Vector{T}
    zs::Vector{T}
    rs::Base.RefValue{T}
    function PauliString{T}(n) where T
        xs = zeros(T, n)
        zs = zeros(T, n)
        rs = Ref(zero(T))
        return new{T}(n, xs, zs, rs)
    end
    function PauliString(n, xs::AbstractVector{T},
                         zs::AbstractVector{T}, rs) where T<:Union{Bool,Base.BitInteger}
        return new{T}(n, xs, zs, Ref(rs))
    end
end
Base.eltype(::Type{PauliString{T}}) where T = T

function Base.empty!(str::PauliString{T}) where T
    str.xs .= zero(T)
    str.zs .= zero(T)
    str.rs[] = zero(T)
    return str
end

Base.@propagate_inbounds @inline function apply!(str::PauliString, gate::Clifford1Q, a)
    @boundscheck check_qubit_bound(str.n, a)
    @inbounds apply!(gate, @view(str.xs[a]), @view(str.zs[a]), str.rs)
    return str
end

Base.@propagate_inbounds @inline function apply!(str::PauliString, gate::Clifford2Q,
                                                 a, b)
    @boundscheck check_qubit_bound(str.n, a)
    @boundscheck check_qubit_bound(str.n, b)
    @inbounds apply!(gate, @view(str.xs[a]), @view(str.zs[a]),
                     @view(str.xs[b]), @view(str.zs[b]), str.rs)
    return str
end

function Base.show(io::IO, str::PauliString{Bool})
    write(io, str.rs[] ? '-' : '+')
    for (x, z) in zip(str.xs, str.zs)
        write(io, x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I'))
    end
end
function Base.show(io::IO, str::PauliString{T}) where T
    for i in 0:(sizeof(T) * 8 - 1)
        write(io, (str.rs[] >> i) != 0 ? '-' : '+')
        for (x, z) in zip(str.xs, str.zs)
            x = (x >> i) != 0
            z = (z >> i) != 0
            write(io, x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I'))
        end
        println(io)
    end
end

# Ignore signs
@inline function inject_pauli!(str::PauliString, x, z, i)
    @boundscheck check_qubit_bound(str.n, i)
    @inbounds str.xs[i] ⊻= x
    @inbounds str.zs[i] ⊻= z
    return str
end

@inline function _cast_bits(::Type{T}, v) where T
    if v isa T
        return v
    elseif v isa Bool
        return v ? (zero(T) - one(T)) : zero(T)
    else
        return v % T
    end
end

function pauli_commute(xas, zas, xbs, zbs)
    TA = eltype(xas)
    TB = eltype(xbs)
    if sizeof(TA) !== sizeof(TB) && TA !== Bool && TB !== Bool
        throw(ArgumentError("Pauli strings dimension mismatch"))
    end
    T = TA === Bool ? TB : TA
    commute = ~zero(T)
    n = length(xas)
    @inbounds @simd for i in 1:n
        xa = _cast_bits(T, xas[i])
        za = _cast_bits(T, zas[i])
        xb = _cast_bits(T, xbs[i])
        zb = _cast_bits(T, zbs[i])
        commute ⊻= (xb & za) ⊻ (xa & zb)
    end
    return commute
end
pauli_commute(stra::PauliString, xbs, zbs) =
    pauli_commute(stra.xs, stra.zs, xbs, zbs)
pauli_commute(xas, zas, strb::PauliString) =
    pauli_commute(xas, zas, strb.xs, strb.zs)
pauli_commute(stra::PauliString, strb::PauliString) =
    pauli_commute(stra.xs, stra.zs, strb.xs, strb.zs)

Base.@propagate_inbounds @inline function pauli_commute_x(stra::PauliString{T},
                                                          xbs) where T
    zas = stra.zs
    commute = ~zero(T)
    for i in xbs
        commute ⊻= zas[i]
    end
    return commute
end
Base.@propagate_inbounds @inline function pauli_commute_z(stra::PauliString{T},
                                                          zbs) where T
    xas = stra.xs
    commute = ~zero(T)
    for i in zbs
        commute ⊻= xas[i]
    end
    return commute
end

# Assume that the expected value without any error is false
function measure_stabilizer(str::PauliString, xs, zs, r=false)
    return ~pauli_commute(str, xs, zs)
end
Base.@propagate_inbounds @inline function measure_stabilizer_x(str::PauliString,
                                                               xs, r=false)
    return ~pauli_commute_x(str, xs)
end
Base.@propagate_inbounds @inline function measure_stabilizer_z(str::PauliString, zs,
                                                               r=false)
    return ~pauli_commute_z(str, zs)
end

const ChT = UInt64
const _chunk_len = sizeof(ChT) * 8
@inline function _get_chunk_mask(bitidx)
    bitidx_0 = bitidx - 1
    chunk = bitidx_0 ÷ _chunk_len + 1
    subbitidx_0 = bitidx_0 % _chunk_len
    mask = one(ChT) << subbitidx_0
    return chunk, mask
end
@inline _getbit(chunk, mask) = chunk & mask != 0
@inline _setbit(chunk, val, mask) = ifelse(val, chunk | mask, chunk & ~mask)

struct StabilizerState
    n::Int
    xs::Matrix{ChT}
    zs::Matrix{ChT}
    rs::Matrix{ChT}
    wxzs::Matrix{ChT}
    function StabilizerState(n)
        nchunks = (2n - 1) ÷ _chunk_len + 1
        xs = zeros(ChT, nchunks, n)
        zs = zeros(ChT, nchunks, n)
        rs = zeros(ChT, nchunks, 1)
        @inbounds for i in 1:n
            chunk1, mask1 = _get_chunk_mask(i)
            xs[chunk1, i] = mask1
            chunk2, mask2 = _get_chunk_mask(i + n)
            zs[chunk2, i] = mask2
        end
        return new(n, xs, zs, rs, Array{ChT}(undef, n, 2))
    end
end
Base.eltype(::Type{StabilizerState}) = Bool

Base.@propagate_inbounds @inline function apply!(state::StabilizerState,
                                                 gate::Clifford1Q, a)
    @boundscheck check_qubit_bound(state.n, a)
    xs = state.xs
    zs = state.zs
    rs = state.rs
    nchunks = size(state.rs, 1)
    @inbounds @simd ivdep for i in 1:nchunks
        apply!(gate, @view(xs[i, a]), @view(zs[i, a]), @view(rs[i, 1]))
    end
    return state
end

Base.@propagate_inbounds @inline function apply!(state::StabilizerState,
                                                 gate::Clifford2Q, a, b)
    @boundscheck check_qubit_bound(state.n, a)
    @boundscheck check_qubit_bound(state.n, b)
    xs = state.xs
    zs = state.zs
    rs = state.rs
    nchunks = size(state.rs, 1)
    @inbounds @simd ivdep for i in 1:nchunks
        apply!(gate, @view(xs[i, a]), @view(zs[i, a]),
               @view(xs[i, b]), @view(zs[i, b]), @view(rs[i, 1]))
    end
    return state
end

Base.@propagate_inbounds @inline function _getindex_x(state::StabilizerState, i, j)
    chunk, mask = _get_chunk_mask(i)
    return _getbit(state.xs[chunk, j], mask)
end

Base.@propagate_inbounds @inline function _getindex_z(state::StabilizerState, i, j)
    chunk, mask = _get_chunk_mask(i)
    return _getbit(state.zs[chunk, j], mask)
end

Base.@propagate_inbounds @inline function _getindex_r(state::StabilizerState, i)
    chunk, mask = _get_chunk_mask(i)
    return _getbit(state.rs[chunk, 1], mask)
end

function Base.show(io::IO, state::StabilizerState)
    for i in 1:state.n
        show(io, get_stabilizer(state, i))
        println(io)
    end
end

function get_stabilizer(state::StabilizerState, i)
    n = state.n
    return PauliString(n, [_getindex_x(state, i + n, j) for j in 1:n],
                       [_getindex_z(state, i + n, j) for j in 1:n],
                       _getindex_r(state, i + n))
end

function init_state_z!(state::StabilizerState, v::Bool=false)
    n = state.n
    xs = state.xs
    zs = state.zs
    xs .= 0
    zs .= 0
    @inbounds for i in 1:n
        chunk1, mask1 = _get_chunk_mask(i)
        xs[chunk1, i] = mask1
        chunk2, mask2 = _get_chunk_mask(i + n)
        zs[chunk2, i] = mask2
    end
    state.rs .= _cast_bits(ChT, v)
    return state
end

function init_state_x!(state::StabilizerState, v::Bool=false)
    n = state.n
    xs = state.xs
    zs = state.zs
    xs .= 0
    zs .= 0
    @inbounds for i in 1:n
        chunk1, mask1 = _get_chunk_mask(i)
        zs[chunk1, i] = mask1
        chunk2, mask2 = _get_chunk_mask(i + n)
        xs[chunk2, i] = mask2
    end
    state.rs .= _cast_bits(ChT, v)
    return state
end

@inline function _accum_pauli_prod_phase(hi, lo, x1, z1, x2, z2)
    v1 = x1 & z2
    v2 = x2 & z1
    m = (z2 ⊻ x1) | ~(x2 | z1)
    change = v1 ⊻ v2
    hi = hi ⊻ ((m ⊻ lo) & change)
    lo = lo ⊻ change
    return hi, lo
end

@inline function _clifford_rowsum1!(state::StabilizerState, chunk_h, mask_h,
                                    chunk_i, mask_i)
    xs = state.xs
    zs = state.zs
    hi = zero(ChT)
    lo = zero(ChT)
    @inbounds for j in 1:state.n
        xi = _getbit(xs[chunk_i, j], mask_i)
        zi = _getbit(zs[chunk_i, j], mask_i)
        if xi
            xh = xs[chunk_h, j]
            zh = zs[chunk_h, j]
            xs[chunk_h, j] = xh ⊻ mask_h
            if zi
                change = zh ⊻ xh
                hi = hi ⊻ ((xh ⊻ lo) & change)
                lo = lo ⊻ change
                zs[chunk_h, j] = zh ⊻ mask_h
            else
                hi = hi ⊻ ((~xh ⊻ lo) & zh)
                lo = lo ⊻ zh
            end
        elseif zi
            xh = xs[chunk_h, j]
            zh = zs[chunk_h, j]
            hi = hi ⊻ ((zh ⊻ lo) & xh)
            lo = lo ⊻ xh
            zs[chunk_h, j] = zh ⊻ mask_h
        end
    end
    @inbounds begin
        rs = state.rs
        rh = rs[chunk_h]
        hi = ifelse(_getbit(rs[chunk_i], mask_i), ~hi, hi)
        rs[chunk_h] = rh ⊻ (hi & mask_h)
    end
    return state
end

@inline function _clifford_rowsum1_2!(state::StabilizerState, chunk_h, mask_h,
                                      chunk_i, mask_i)
    xs = state.xs
    zs = state.zs
    hi = zero(ChT)
    lo = zero(ChT)
    @inbounds for j in 1:state.n
        xi = _getbit(xs[chunk_i, j], mask_i)
        zi = _getbit(zs[chunk_i, j], mask_i)
        xh = xs[chunk_h, j]
        zh = zs[chunk_h, j]

        v1 = ifelse(xi, zh, zero(ChT))
        v2 = ifelse(zi, xh, zero(ChT))
        m = ifelse(xi, ~zh, zh) | ifelse(zi, zero(ChT), ~xh)
        change = v1 ⊻ v2
        hi = hi ⊻ ((m ⊻ lo) & change)
        lo = lo ⊻ change

        xs[chunk_h, j] = ifelse(xi, xh ⊻ mask_h, xh)
        zs[chunk_h, j] = ifelse(zi, zh ⊻ mask_h, zh)
    end
    @inbounds begin
        rs = state.rs
        rh = rs[chunk_h]
        hi = ifelse(_getbit(rs[chunk_i], mask_i), ~hi, hi)
        rs[chunk_h] = rh ⊻ (hi & mask_h)
    end
    return state
end
@inline function _clifford_rowsum2!(state::StabilizerState, chunk_i, mask_i, wrs)
    xs = state.xs
    zs = state.zs
    wxzs = state.wxzs
    hi = zero(ChT)
    lo = zero(ChT)
    @inbounds for j in 1:state.n
        xi = xs[chunk_i, j]
        xh = wxzs[j, 1]
        zi = zs[chunk_i, j]
        zh = wxzs[j, 2]
        hi, lo = _accum_pauli_prod_phase(hi, lo, xi, zi, xh, zh)
        wxzs[j, 1] = xh ⊻ (xi & mask_i)
        wxzs[j, 2] = zh ⊻ (zi & mask_i)
    end
    @inbounds begin
        wrs[] ⊻= (hi ⊻ state.rs[chunk_i, 1]) & mask_i
    end
    return state
end
@inline function _clifford_rowsum3!(state::StabilizerState, wrs)
    @assert sizeof(ChT) == 8
    wxzs = state.wxzs
    hi = zero(ChT)
    lo = zero(ChT)
    @inbounds for j in 1:state.n
        x = wxzs[j, 1]
        z = wxzs[j, 2]

        x1 = (x >> 32) % UInt32
        z1 = (z >> 32) % UInt32
        x2 = x % UInt32
        z2 = z % UInt32

        xh = x1 % UInt64
        zh = z1 % UInt64
        xi = x2 % UInt64
        zi = z2 % UInt64

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        x1 = x >> 16
        z1 = z >> 16
        x2 = x & 0x0000ffff
        z2 = z & 0x0000ffff

        xh |= (x1 % UInt64) << 32
        zh |= (z1 % UInt64) << 32
        xi |= (x2 % UInt64) << 32
        zi |= (z2 % UInt64) << 32

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        x1 = x >> 8
        z1 = z >> 8
        x2 = x & 0x000000ff
        z2 = z & 0x000000ff

        xh |= (x1 % UInt64) << 48
        zh |= (z1 % UInt64) << 48
        xi |= (x2 % UInt64) << 48
        zi |= (z2 % UInt64) << 48

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        x1 = x >> 4
        z1 = z >> 4
        x2 = x & 0x0000000f
        z2 = z & 0x0000000f

        xh |= (x1 % UInt64) << 56
        zh |= (z1 % UInt64) << 56
        xi |= (x2 % UInt64) << 56
        zi |= (z2 % UInt64) << 56

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        x1 = x >> 2
        z1 = z >> 2
        x2 = x & 0x00000003
        z2 = z & 0x00000003

        xh |= (x1 % UInt64) << 60
        zh |= (z1 % UInt64) << 60
        xi |= (x2 % UInt64) << 60
        zi |= (z2 % UInt64) << 60

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        x1 = x >> 1
        z1 = z >> 1
        x2 = x & 0x00000001
        z2 = z & 0x00000001

        xh |= (x1 % UInt64) << 62
        zh |= (z1 % UInt64) << 62
        xi |= (x2 % UInt64) << 62
        zi |= (z2 % UInt64) << 62

        hi, lo = _accum_pauli_prod_phase(hi, lo, xh, zh, xi, zi)
    end
    return (count_ones(wrs[] ⊻ hi) ⊻ (count_ones(lo) >> 1)) & 1 != 0
end

@inline function _inject_zresult!(state::StabilizerState, n, a, p, chunk2, mask2, res)
    chunk1, mask1 = _get_chunk_mask(p - n)
    xs = state.xs
    zs = state.zs
    @inbounds for i in 1:n
        # state.xs[i][p - n] = state.xs[i][p]
        # state.xs[i][p] = false
        x2 = xs[chunk2, i]
        xs[chunk2, i] = x2 & ~mask2
        # Note that the load of x1 has to happen after the store of x2
        # since the two may alias
        xs[chunk1, i] = _setbit(xs[chunk1, i], (x2 & mask2) != 0, mask1)

        # state.zs[i][p - n] = state.zs[i][p]
        # state.zs[i][p] = i == a
        z2 = zs[chunk2, i]
        zs[chunk2, i] = _setbit(z2, i == a, mask2)
        # Note that the load of z1 has to happen after the store of z2
        # since the two may alias
        zs[chunk1, i] = _setbit(zs[chunk1, i], (z2 & mask2) != 0, mask1)
    end
    # state.rs[p - n] = state.rs[p]
    # state.rs[p] = res
    rs = state.rs
    @inbounds begin
        r2 = rs[chunk2, 1]
        rs[chunk2, 1] = _setbit(r2, res, mask2)
        # Note that the load of r1 has to happen after the store of r2
        # since the two may alias
        rs[chunk1, 1] = _setbit(rs[chunk1, 1], (r2 & mask2) != 0, mask1)
    end
    return
end

# Randomly pick a result
function measure_z!(state::StabilizerState, a; force=nothing)
    n = state.n
    check_qubit_bound(n, a)
    p = 0
    chunk_p = 0
    mask_p = zero(ChT)
    found_p = false
    @inbounds for _p in (n + 1):(2 * n)
        chunk_p, mask_p = _get_chunk_mask(_p)
        if _getbit(state.xs[chunk_p, a], mask_p)
            p = _p
            found_p = true
            break
        end
    end
    @inbounds if found_p
        nchunks = size(state.rs, 1)
        if nchunks <= 1
            for i in 1:nchunks
                mask_i = state.xs[i, a]
                bitidx = (p - 1 - (i - 1) * _chunk_len) % UInt
                mask_i = ifelse(bitidx < _chunk_len,
                                mask_i & ~(one(ChT) << bitidx), mask_i)
                if mask_i == 0
                    continue
                end
                _clifford_rowsum1_2!(state, i, mask_i, chunk_p, mask_p)
            end
        else
            for i in 1:nchunks
                mask_i = state.xs[i, a]
                bitidx = (p - 1 - (i - 1) * _chunk_len) % UInt
                mask_i = ifelse(bitidx < _chunk_len,
                                mask_i & ~(one(ChT) << bitidx), mask_i)
                if mask_i == 0
                    continue
                end
                _clifford_rowsum1!(state, i, mask_i, chunk_p, mask_p)
            end
        end
        res = force !== nothing ? force : rand(Bool)
        _inject_zresult!(state, n, a, p, chunk_p, mask_p, res)
        return res, false
    else
        nfullchunks, rembits = divrem(n, _chunk_len)
        state.wxzs .= zero(ChT)
        wrs = Ref(zero(ChT))
        total_mask = zero(ChT)
        if rembits == 0
            for i in 1:nfullchunks
                mask_i = state.xs[i, a]
                if mask_i != 0
                    total_mask |= mask_i
                    _clifford_rowsum2!(state, i + nfullchunks, mask_i, wrs)
                end
            end
        else
            mask = zero(ChT)
            for i in 1:nfullchunks
                new_mask = state.xs[i, a]
                mask_i = mask | (new_mask << rembits)
                mask = new_mask >> (_chunk_len - rembits)
                if mask_i != 0
                    total_mask |= mask_i
                    _clifford_rowsum2!(state, i + nfullchunks, mask_i, wrs)
                end
            end
            new_mask = state.xs[nfullchunks + 1, a]
            mask_i = mask | (new_mask << rembits)
            rembits2 = rembits * 2
            mask_i &= (one(ChT) << rembits2) - one(ChT)
            if mask_i != 0
                total_mask |= mask_i
                _clifford_rowsum2!(state, 2 * nfullchunks + 1, mask_i, wrs)
            end
            if rembits2 > _chunk_len
                mask_i = new_mask >> (_chunk_len - rembits)
                # Mask out the end of the range
                mask_i &= (one(ChT) << (rembits2 - _chunk_len)) - one(ChT)
                if mask_i != 0
                    total_mask |= mask_i
                    _clifford_rowsum2!(state, 2 * nfullchunks + 2, mask_i, wrs)
                end
            end
        end
        if count_ones(total_mask) == 1
            return wrs[] != 0, true
        end
        return _clifford_rowsum3!(state, wrs), true
    end
end

struct InvStabilizerState
    n::Int
    xzs::Array{ChT,3}
    rs::Matrix{Bool}
    function InvStabilizerState(n)
        # Align the chunk number to 2
        nchunks = ((n - 1) ÷ (_chunk_len * 2) + 1) * 2
        # XX,            XZ,            ZX,            ZZ
        # X stab X term, X stab Z term, Z stab X term, Z stab Z term
        xzs = zeros(ChT, nchunks, 4, n)
        rs = zeros(Bool, n, 2)
        @inbounds for i in 1:n
            chunk, mask = _get_chunk_mask(i)
            xzs[chunk, 1, i] = mask
            xzs[chunk, 4, i] = mask
        end
        return new(n, xzs, rs)
    end
end
Base.eltype(::Type{InvStabilizerState}) = Bool

function Base.show(io::IO, state::InvStabilizerState)
    for i in 1:state.n
        print(io, "X[$i]: ")
        show(io, get_inv_stabilizer(state, i, false))
        println(io)
        print(io, "Z[$i]: ")
        show(io, get_inv_stabilizer(state, i, true))
        println(io)
    end
end

Base.@propagate_inbounds @inline function _getindex(state::InvStabilizerState,
                                                    i, j, k)
    chunk, mask = _get_chunk_mask(i)
    return _getbit(state.xzs[chunk, j, k], mask)
end

function get_inv_stabilizer(state::InvStabilizerState, i, z::Bool)
    n = state.n
    k = z ? 2 : 1
    return PauliString(n, [_getindex(state, j, 2k - 1, i) for j in 1:n],
                       [_getindex(state, j, 2k, i) for j in 1:n],
                       state.rs[i, k])
end

function measure_x!(state::StabilizerState, a; force=nothing)
    @boundscheck check_qubit_bound(state.n, a)
    @inbounds apply!(state, HGate(), a)
    res = measure_z!(state, a; force=force)
    @inbounds apply!(state, HGate(), a)
    return res
end

function measure_y!(state::StabilizerState, a; force=nothing)
    @boundscheck check_qubit_bound(state.n, a)
    @inbounds apply!(state, SXGate(), a)
    res = measure_z!(state, a; force=force)
    @inbounds apply!(state, ISXGate(), a)
    return res
end

function measure_zs!(state::StabilizerState, idxs; force=nothing)
    if isempty(idxs)
        return false, true
    end
    idx0 = idxs[1]
    nidxs = length(idxs)
    if nidxs == 1
        return measure_z!(state, idx0; force=force)
    end
    @boundscheck check_qubit_bound(state.n, idx0)
    for i in 2:nidxs
        idxi = @inbounds idxs[i]
        @boundscheck check_qubit_bound(state.n, idxi)
        apply!(state, CNOTGate(), idxi, idx0)
    end
    res = measure_z!(state, idx0; force=force)
    @inbounds for i in 2:nidxs
        apply!(state, CNOTGate(), idxs[i], idx0)
    end
    return res
end

function measure_xs!(state::StabilizerState, idxs; force=nothing)
    for idx in idxs
        apply!(state, HGate(), idx)
    end
    res = @inbounds measure_zs!(state, idxs; force=force)
    @inbounds for idx in idxs
        apply!(state, HGate(), idx)
    end
    return res
end

function measure_ys!(state::StabilizerState, idxs; force=nothing)
    for idx in idxs
        apply!(state, SXGate(), idx)
    end
    res = @inbounds measure_zs!(state, idxs; force=force)
    @inbounds for idx in idxs
        apply!(state, ISXGate(), idx)
    end
    return res
end

function measure_paulis!(state::StabilizerState, xs, zs; force=nothing)
    @assert length(xs) == state.n
    @assert length(zs) == state.n
    local idx0
    for i in 1:state.n
        x = xs[i]
        z = zs[i]
        @inbounds if x
            if z
                apply!(state, SXGate(), i)
            else
                apply!(state, HGate(), i)
            end
        elseif !z
            continue
        end
        if !@isdefined(idx0)
            idx0 = i
        else
            @inbounds apply!(state, CNOTGate(), i, idx0)
        end
    end
    if !@isdefined(idx0)
        return false, true
    end
    res = measure_z!(state, idx0; force=force)
    @inbounds for i in state.n:-1:1
        x = xs[i]
        z = zs[i]
        if !x && !z
            continue
        end
        if i != idx0
            apply!(state, CNOTGate(), i, idx0)
        end
        if x
            if z
                apply!(state, ISXGate(), i)
            else
                apply!(state, HGate(), i)
            end
        end
    end
    return res
end

function measure_paulis!(state::StabilizerState, str::PauliString{Bool}; force=nothing)
    if force !== nothing
        force ⊻= str.rs[]
    end
    v, det = measure_paulis!(state, str.xs, str.ys; force=force)
    v ⊻= str.rs[]
    return v, det
end

function measure_stabilizer(state::StabilizerState, xs, zs, r=false)
    v, det = measure_paulis!(state, xs, ys)
    return v ⊻ r
end
function measure_stabilizer_x(state::StabilizerState, xs, r=false)
    v, det = measure_xs!(state, xs)
    return v ⊻ r
end
function measure_stabilizer_z(state::StabilizerState, zs, r=false)
    v, det = measure_zs!(state, zs)
    return v ⊻ r
end

Base.@propagate_inbounds @inline function inject_pauli!(state::StabilizerState,
                                                        x::Bool, z::Bool, i)
    if x
        if z
            apply!(state, YGate(), i)
        else
            apply!(state, XGate(), i)
        end
    elseif z
        apply!(state, ZGate(), i)
    end
    return state
end

end

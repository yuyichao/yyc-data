#!/usr/bin/julia

module Clifford

using Random
using ..Utils
using SIMD

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

@inline function _cast_bits(::Type{T}, v) where T
    if v isa T
        return v
    elseif v isa Bool
        return v ? (zero(T) - one(T)) : zero(T)
    else
        return v % T
    end
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

@inline vcount_ones_u8(v::Vec{N,T}) where {N,T} =
    count_ones(reinterpret(Vec{N * sizeof(T),UInt8}, v))

@inline function assume(v::Bool)
    Base.llvmcall(
        ("""
         declare void @llvm.assume(i1)
         define void @fw_assume(i8) alwaysinline
         {
             %v = trunc i8 %0 to i1
             call void @llvm.assume(i1 %v)
             ret void
         }
         """, "fw_assume"), Cvoid, Tuple{Bool}, v)
end

@generated function _gep_array(ptr::Ptr{T}, szs, index::NTuple{N}) where {T,N}
    exp = :(index[$N] - 1)
    for i in (N - 1):-1:1
        exp = :($exp * szs[$i] + index[$i] - 1)
    end
    return :(@inline; ptr + $exp * $(sizeof(T)))
end

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
        write(io, (str.rs[] >> i) & 1 != 0 ? '-' : '+')
        for (x, z) in zip(str.xs, str.zs)
            x = (x >> i) & 1 != 0
            z = (z >> i) & 1 != 0
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
    assume(nchunks > 0)
    @inbounds @simd ivdep for i in 1:nchunks
        apply!(gate, @view(xs[i, a]), @view(zs[i, a]), @view(rs[i, 1]))
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::StabilizerState,
                                                 gate::IGate, a)
    @boundscheck check_qubit_bound(state.n, a)
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
    assume(nchunks > 0)
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

function init_state_z!(state::StabilizerState, v::Bool)
    n = state.n
    xs = state.xs
    zs = state.zs
    xs .= 0
    zs .= 0
    assume(n > 0)
    @inbounds for i in 1:n
        chunk1, mask1 = _get_chunk_mask(i)
        xs[chunk1, i] = mask1
        chunk2, mask2 = _get_chunk_mask(i + n)
        zs[chunk2, i] = mask2
    end
    state.rs .= _cast_bits(ChT, v)
    return state
end

function init_state_x!(state::StabilizerState, v::Bool)
    n = state.n
    xs = state.xs
    zs = state.zs
    xs .= 0
    zs .= 0
    assume(n > 0)
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
    assume(state.n > 0)
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
    assume(state.n > 0)
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
    assume(state.n > 0)
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
    assume(state.n > 0)
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
    assume(n > 0)
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

Base.@propagate_inbounds @inline function measure_z!(state::StabilizerState, a;
                                                     force=nothing)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    return _measure_z!(state, n, a, force)
end

# Randomly pick a result
function _measure_z!(state::StabilizerState, n, a, force)
    assume(n == state.n)
    p = 0
    chunk_p = 0
    mask_p = zero(ChT)
    found_p = false
    assume(n > 0)
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
        assume(nchunks > 0)
        if nchunks <= 1
            mask_i = state.xs[1, a]
            bitidx = (p - 1) % UInt
            mask_i = ifelse(bitidx < _chunk_len,
                            mask_i & ~(one(ChT) << bitidx), mask_i)
            if mask_i != 0
                _clifford_rowsum1_2!(state, 1, mask_i, chunk_p, mask_p)
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
            assume(nfullchunks > 0)
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
    ws::Matrix{ChT}
    function InvStabilizerState(n)
        # Align the chunk number to 2
        nchunks = ((n - 1) ÷ (_chunk_len * 2) + 1) * 2
        # XX,            XZ,            ZX,            ZZ
        # X stab X term, X stab Z term, Z stab X term, Z stab Z term
        xzs = zeros(ChT, nchunks, 4, n)
        rs = zeros(Bool, n, 2)
        assume(size(xzs, 2) == 4)
        @inbounds for i in 1:n
            chunk, mask = _get_chunk_mask(i)
            xzs[chunk, 1, i] = mask
            xzs[chunk, 4, i] = mask
        end
        return new(n, xzs, rs, zeros(ChT, 4, nchunks >> 1))
        # return new(n, xzs, rs, zeros(ChT, 6 * n, 2))
    end
end
Base.eltype(::Type{InvStabilizerState}) = Bool

function Base.show(io::IO, state::InvStabilizerState)
    assume(state.n > 0)
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
    xzs = state.xzs
    assume(size(xzs, 2) == 4)
    return _getbit(xzs[chunk, j, k], mask)
end

function get_inv_stabilizer(state::InvStabilizerState, i, z::Bool)
    n = state.n
    k = z ? 2 : 1
    assume(n > 0)
    return PauliString(n, [_getindex(state, j, 2k - 1, i) for j in 1:n],
                       [_getindex(state, j, 2k, i) for j in 1:n],
                       state.rs[i, k])
end

function init_state_z!(state::InvStabilizerState, v::Bool)
    n = state.n
    xzs = state.xzs
    assume(size(xzs, 2) == 4)
    xzs .= 0
    assume(n > 0)
    @inbounds for i in 1:n
        chunk, mask = _get_chunk_mask(i)
        xzs[chunk, 1, i] = mask
        xzs[chunk, 4, i] = mask
    end
    rs = state.rs
    assume(size(rs, 1) == n)
    @inbounds begin
        rs[:, 1] .= false
        rs[:, 2] .= v
    end
    return state
end

function init_state_x!(state::InvStabilizerState, v::Bool)
    n = state.n
    xzs = state.xzs
    assume(size(xzs, 2) == 4)
    xzs .= 0
    assume(n > 0)
    @inbounds for i in 1:n
        chunk, mask = _get_chunk_mask(i)
        xzs[chunk, 2, i] = mask
        xzs[chunk, 3, i] = mask
    end
    rs = state.rs
    assume(size(rs, 1) == n)
    @inbounds begin
        rs[:, 1] .= v
        rs[:, 2] .= false
    end
    return state
end

# P1 = P1 * P2
@inline function pauli_multiply!(px1s, pz1s, px2s, pz2s, n)
    VT8 = Vec{8,ChT}
    VT4 = Vec{4,ChT}
    VT2 = Vec{2,ChT}

    clo = zero(VT2)
    chi = zero(VT2)
    # Manual jump threading for the small n cases
    if n <= 2
        xz1 = vload(VT4, px1s)
        xz2 = vload(VT4, px2s)
        new_xz1 = xz1 ⊻ xz2
        vstore(new_xz1, px1s)

        x1 = Vec((xz1[1], xz1[2]))
        z1 = Vec((xz1[3], xz1[4]))
        x2 = Vec((xz2[1], xz2[2]))
        z2 = Vec((xz2[3], xz2[4]))
        new_x1 = Vec((new_xz1[1], new_xz1[2]))
        new_z1 = Vec((new_xz1[3], new_xz1[4]))

        v1 = x1 & z2
        v2 = x2 & z1
        m = new_x1 ⊻ new_z1 ⊻ v1
        change = v1 ⊻ v2
        chi = m & change
        clo = change
        @goto n2_case_end
    elseif n < 8
        @goto n4_case
    end
    hi = zero(VT8)
    lo = zero(VT8)
    nalign = n & ~7
    # LLVM may think the vstore in the loop aliases the array pointer
    # so we extract the array pointer out of the loop.
    @inbounds for i0 in 1:(n >> 3)
        i = (i0 - 1) * 8 + 1
        x1 = vload(VT8, _gep_array(px1s, (), (i,)))
        x2 = vload(VT8, _gep_array(px2s, (), (i,)))
        new_x1 = x1 ⊻ x2
        vstore(new_x1, _gep_array(px1s, (), (i,)))
        z1 = vload(VT8, _gep_array(pz1s, (), (i,)))
        z2 = vload(VT8, _gep_array(pz2s, (), (i,)))
        new_z1 = z1 ⊻ z2
        vstore(new_z1, _gep_array(pz1s, (), (i,)))

        v1 = x1 & z2
        v2 = x2 & z1
        m = new_x1 ⊻ new_z1 ⊻ v1
        change = v1 ⊻ v2
        hi = hi ⊻ ((m ⊻ lo) & change)
        lo = lo ⊻ change
    end
    clo = VT2((lo[1] ⊻ lo[3] ⊻ lo[5] ⊻ lo[7],
               lo[2] ⊻ lo[4] ⊻ lo[6] ⊻ lo[8]))
    chi = VT2((hi[1] ⊻ hi[3] ⊻ hi[5] ⊻ hi[7],
               hi[2] ⊻ hi[4] ⊻ hi[6] ⊻ hi[8]))
    px1s += nalign * 8
    pz1s += nalign * 8
    px2s += nalign * 8
    pz2s += nalign * 8

    if n & 4 == 0
        @goto n4_case_end
    end
    @label n4_case
    x1 = vload(VT4, px1s)
    x2 = vload(VT4, px2s)
    new_x1 = x1 ⊻ x2
    vstore(new_x1, px1s)
    z1 = vload(VT4, pz1s)
    z2 = vload(VT4, pz2s)
    new_z1 = z1 ⊻ z2
    vstore(new_z1, pz1s)

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    hi4 = m & change
    lo4 = change
    clo ⊻= VT2((lo4[1] ⊻ lo4[3], lo4[2] ⊻ lo4[4]))
    chi ⊻= VT2((hi4[1] ⊻ hi4[3], hi4[2] ⊻ hi4[4]))
    px1s += sizeof(VT4)
    pz1s += sizeof(VT4)
    px2s += sizeof(VT4)
    pz2s += sizeof(VT4)
    @label n4_case_end

    if n & 2 == 0
        @goto n2_case_end
    end
    @label n2_case
    x1 = vloada(VT2, px1s)
    x2 = vloada(VT2, px2s)
    new_x1 = x1 ⊻ x2
    vstorea(new_x1, px1s)
    z1 = vloada(VT2, pz1s)
    z2 = vloada(VT2, pz2s)
    new_z1 = z1 ⊻ z2
    vstorea(new_z1, pz1s)

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    chi ⊻= m & change
    clo ⊻= change
    @label n2_case_end

    cnt = vcount_ones_u8(clo) + vcount_ones_u8(chi) << 1
    return reduce(+, cnt) & 0x3
end

@inline function pauli_multiply_2!(px1s_1, pz1s_1, px2s_1, pz2s_1,
                                   px1s_2, pz1s_2, px2s_2, pz2s_2, n)
    VT8 = Vec{8,ChT}
    VT4 = Vec{4,ChT}
    VT2 = Vec{2,ChT}

    clo_1 = zero(VT2)
    chi_1 = zero(VT2)
    clo_2 = zero(VT2)
    chi_2 = zero(VT2)
    # Manual jump threading for the small n cases
    if n <= 2
        xz1 = vload(VT4, px1s_1)
        xz2 = vload(VT4, px2s_1)
        new_xz1 = xz1 ⊻ xz2
        vstore(new_xz1, px1s_1)

        x1 = Vec((xz1[1], xz1[2]))
        z1 = Vec((xz1[3], xz1[4]))
        x2 = Vec((xz2[1], xz2[2]))
        z2 = Vec((xz2[3], xz2[4]))
        new_x1 = Vec((new_xz1[1], new_xz1[2]))
        new_z1 = Vec((new_xz1[3], new_xz1[4]))

        v1 = x1 & z2
        v2 = x2 & z1
        m = new_x1 ⊻ new_z1 ⊻ v1
        change = v1 ⊻ v2
        chi_1 = m & change
        clo_1 = change

        xz1 = vload(VT4, px1s_2)
        xz2 = vload(VT4, px2s_2)
        new_xz1 = xz1 ⊻ xz2
        vstore(new_xz1, px1s_2)

        x1 = Vec((xz1[1], xz1[2]))
        z1 = Vec((xz1[3], xz1[4]))
        x2 = Vec((xz2[1], xz2[2]))
        z2 = Vec((xz2[3], xz2[4]))
        new_x1 = Vec((new_xz1[1], new_xz1[2]))
        new_z1 = Vec((new_xz1[3], new_xz1[4]))

        v1 = x1 & z2
        v2 = x2 & z1
        m = new_x1 ⊻ new_z1 ⊻ v1
        change = v1 ⊻ v2
        chi_2 = m & change
        clo_2 = change
        @goto n2_case_end
    elseif n < 8
        @goto n4_case
    end
    hi_1 = zero(VT8)
    lo_1 = zero(VT8)
    hi_2 = zero(VT8)
    lo_2 = zero(VT8)
    nalign = n & ~7
    # LLVM may think the vstore in the loop aliases the array pointer
    # so we extract the array pointer out of the loop.
    @inbounds for i0 in 1:(n >> 3)
        i = (i0 - 1) * 8 + 1
        x1 = vload(VT8, _gep_array(px1s_1, (), (i,)))
        x2 = vload(VT8, _gep_array(px2s_1, (), (i,)))
        new_x1 = x1 ⊻ x2
        vstore(new_x1, _gep_array(px1s_1, (), (i,)))
        z1 = vload(VT8, _gep_array(pz1s_1, (), (i,)))
        z2 = vload(VT8, _gep_array(pz2s_1, (), (i,)))
        new_z1 = z1 ⊻ z2
        vstore(new_z1, _gep_array(pz1s_1, (), (i,)))

        v1 = x1 & z2
        v2 = x2 & z1
        m = new_x1 ⊻ new_z1 ⊻ v1
        change = v1 ⊻ v2
        hi_1 = hi_1 ⊻ ((m ⊻ lo_1) & change)
        lo_1 = lo_1 ⊻ change

        x1 = vload(VT8, _gep_array(px1s_2, (), (i,)))
        x2 = vload(VT8, _gep_array(px2s_2, (), (i,)))
        new_x1 = x1 ⊻ x2
        vstore(new_x1, _gep_array(px1s_2, (), (i,)))
        z1 = vload(VT8, _gep_array(pz1s_2, (), (i,)))
        z2 = vload(VT8, _gep_array(pz2s_2, (), (i,)))
        new_z1 = z1 ⊻ z2
        vstore(new_z1, _gep_array(pz1s_2, (), (i,)))

        v1 = x1 & z2
        v2 = x2 & z1
        m = new_x1 ⊻ new_z1 ⊻ v1
        change = v1 ⊻ v2
        hi_2 = hi_2 ⊻ ((m ⊻ lo_2) & change)
        lo_2 = lo_2 ⊻ change
    end
    clo_1 = VT2((lo_1[1] ⊻ lo_1[3] ⊻ lo_1[5] ⊻ lo_1[7],
                 lo_1[2] ⊻ lo_1[4] ⊻ lo_1[6] ⊻ lo_1[8]))
    chi_1 = VT2((hi_1[1] ⊻ hi_1[3] ⊻ hi_1[5] ⊻ hi_1[7],
                 hi_1[2] ⊻ hi_1[4] ⊻ hi_1[6] ⊻ hi_1[8]))
    clo_2 = VT2((lo_2[1] ⊻ lo_2[3] ⊻ lo_2[5] ⊻ lo_2[7],
                 lo_2[2] ⊻ lo_2[4] ⊻ lo_2[6] ⊻ lo_2[8]))
    chi_2 = VT2((hi_2[1] ⊻ hi_2[3] ⊻ hi_2[5] ⊻ hi_2[7],
                 hi_2[2] ⊻ hi_2[4] ⊻ hi_2[6] ⊻ hi_2[8]))
    px1s_1 += nalign * 8
    pz1s_1 += nalign * 8
    px2s_1 += nalign * 8
    pz2s_1 += nalign * 8
    px1s_2 += nalign * 8
    pz1s_2 += nalign * 8
    px2s_2 += nalign * 8
    pz2s_2 += nalign * 8

    if n & 4 == 0
        @goto n4_case_end
    end
    @label n4_case
    x1 = vload(VT4, px1s_1)
    x2 = vload(VT4, px2s_1)
    new_x1 = x1 ⊻ x2
    vstore(new_x1, px1s_1)
    z1 = vload(VT4, pz1s_1)
    z2 = vload(VT4, pz2s_1)
    new_z1 = z1 ⊻ z2
    vstore(new_z1, pz1s_1)

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    hi4 = m & change
    lo4 = change
    clo_1 ⊻= VT2((lo4[1] ⊻ lo4[3], lo4[2] ⊻ lo4[4]))
    chi_1 ⊻= VT2((hi4[1] ⊻ hi4[3], hi4[2] ⊻ hi4[4]))

    px1s_1 += sizeof(VT4)
    pz1s_1 += sizeof(VT4)
    px2s_1 += sizeof(VT4)
    pz2s_1 += sizeof(VT4)

    x1 = vload(VT4, px1s_2)
    x2 = vload(VT4, px2s_2)
    new_x1 = x1 ⊻ x2
    vstore(new_x1, px1s_2)
    z1 = vload(VT4, pz1s_2)
    z2 = vload(VT4, pz2s_2)
    new_z1 = z1 ⊻ z2
    vstore(new_z1, pz1s_2)

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    hi4 = m & change
    lo4 = change
    clo_2 ⊻= VT2((lo4[1] ⊻ lo4[3], lo4[2] ⊻ lo4[4]))
    chi_2 ⊻= VT2((hi4[1] ⊻ hi4[3], hi4[2] ⊻ hi4[4]))

    px1s_2 += sizeof(VT4)
    pz1s_2 += sizeof(VT4)
    px2s_2 += sizeof(VT4)
    pz2s_2 += sizeof(VT4)
    @label n4_case_end

    if n & 2 == 0
        @goto n2_case_end
    end
    @label n2_case
    x1 = vloada(VT2, px1s_1)
    x2 = vloada(VT2, px2s_1)
    new_x1 = x1 ⊻ x2
    vstorea(new_x1, px1s_1)
    z1 = vloada(VT2, pz1s_1)
    z2 = vloada(VT2, pz2s_1)
    new_z1 = z1 ⊻ z2
    vstorea(new_z1, pz1s_1)

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    chi_1 ⊻= m & change
    clo_1 ⊻= change

    x1 = vloada(VT2, px1s_2)
    x2 = vloada(VT2, px2s_2)
    new_x1 = x1 ⊻ x2
    vstorea(new_x1, px1s_2)
    z1 = vloada(VT2, pz1s_2)
    z2 = vloada(VT2, pz2s_2)
    new_z1 = z1 ⊻ z2
    vstorea(new_z1, pz1s_2)

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    chi_2 ⊻= m & change
    clo_2 ⊻= change
    @label n2_case_end

    cnt_1 = vcount_ones_u8(clo_1) + vcount_ones_u8(chi_1) << 1
    cnt_2 = vcount_ones_u8(clo_2) + vcount_ones_u8(chi_2) << 1
    return reduce(+, cnt_1) & 0x3, reduce(+, cnt_2) & 0x3
end

Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::Composite1Q, a)
    apply!(state, gate.g1, a)
    apply!(state, gate.g2, a)
    return state
end
@inline function apply!(state::InvStabilizerState, gate::IGate, a)
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::HGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)
    @inbounds begin
        @simd ivdep for i in 1:nchunks
            xx = xzs[i, 1, a]
            xz = xzs[i, 2, a]
            zx = xzs[i, 3, a]
            zz = xzs[i, 4, a]
            xzs[i, 1, a] = zx
            xzs[i, 2, a] = zz
            xzs[i, 3, a] = xx
            xzs[i, 4, a] = xz
        end
        rx = rs[a, 1]
        rz = rs[a, 2]
        rs[a, 2] = rx
        rs[a, 1] = rz
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::XGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    rs = state.rs
    assume(size(rs, 1) == n)
    @inbounds begin
        rs[a, 2] = ~rs[a, 2]
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::YGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    rs = state.rs
    assume(size(rs, 1) == n)
    @inbounds begin
        rs[a, 1] = ~rs[a, 1]
        rs[a, 2] = ~rs[a, 2]
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::ZGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    rs = state.rs
    assume(size(rs, 1) == n)
    @inbounds begin
        rs[a, 1] = ~rs[a, 1]
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::SGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)
    @inbounds GC.@preserve xzs begin
        px1s = pointer(@view(xzs[1, 1, a]))
        pz1s = pointer(@view(xzs[1, 2, a]))
        px2s = pointer(@view(xzs[1, 3, a]))
        pz2s = pointer(@view(xzs[1, 4, a]))
        prod_phase = pauli_multiply!(px1s, pz1s, px2s, pz2s, nchunks)
        assume(prod_phase & 0x1 != 0)
        rs[a, 1] ⊻= rs[a, 2] ⊻ ((prod_phase & 0x2) != 0)
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::ISGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)
    @inbounds GC.@preserve xzs begin
        px1s = pointer(@view(xzs[1, 1, a]))
        pz1s = pointer(@view(xzs[1, 2, a]))
        px2s = pointer(@view(xzs[1, 3, a]))
        pz2s = pointer(@view(xzs[1, 4, a]))
        prod_phase = pauli_multiply!(px1s, pz1s, px2s, pz2s, nchunks)
        assume(prod_phase & 0x1 != 0)
        rs[a, 1] ⊻= rs[a, 2] ⊻ ((prod_phase & 0x2) == 0)
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::SXGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)
    @inbounds GC.@preserve xzs begin
        px1s = pointer(@view(xzs[1, 3, a]))
        pz1s = pointer(@view(xzs[1, 4, a]))
        px2s = pointer(@view(xzs[1, 1, a]))
        pz2s = pointer(@view(xzs[1, 2, a]))
        prod_phase = pauli_multiply!(px1s, pz1s, px2s, pz2s, nchunks)
        assume(prod_phase & 0x1 != 0)
        rs[a, 2] ⊻= rs[a, 1] ⊻ ((prod_phase & 0x2) != 0)
    end
    return state
end
Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::ISXGate, a)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)
    @inbounds GC.@preserve xzs begin
        px1s = pointer(@view(xzs[1, 3, a]))
        pz1s = pointer(@view(xzs[1, 4, a]))
        px2s = pointer(@view(xzs[1, 1, a]))
        pz2s = pointer(@view(xzs[1, 2, a]))
        prod_phase = pauli_multiply!(px1s, pz1s, px2s, pz2s, nchunks)
        assume(prod_phase & 0x1 != 0)
        rs[a, 2] ⊻= rs[a, 1] ⊻ ((prod_phase & 0x2) == 0)
    end
    return state
end

Base.@propagate_inbounds @inline function apply!(state::InvStabilizerState,
                                                 gate::CNOTGate, a, b)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    @boundscheck check_qubit_bound(n, b)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)
    @inbounds GC.@preserve xzs begin
        # Xa′ = Xa * Xb
        px1s_1 = pointer(@view(xzs[1, 1, a]))
        pz1s_1 = pointer(@view(xzs[1, 2, a]))
        px2s_1 = pointer(@view(xzs[1, 1, b]))
        pz2s_1 = pointer(@view(xzs[1, 2, b]))

        # Zb′ = Za * Zb
        px1s_2 = pointer(@view(xzs[1, 3, b]))
        pz1s_2 = pointer(@view(xzs[1, 4, b]))
        px2s_2 = pointer(@view(xzs[1, 3, a]))
        pz2s_2 = pointer(@view(xzs[1, 4, a]))

        prod_phase_1, prod_phase_2 =
            pauli_multiply_2!(px1s_1, pz1s_1, px2s_1, pz2s_1,
                              px1s_2, pz1s_2, px2s_2, pz2s_2, nchunks)
        assume(prod_phase_1 & 0x1 == 0)
        assume(prod_phase_2 & 0x1 == 0)
        rs[a, 1] ⊻= rs[b, 1] ⊻ (prod_phase_1 != 0)
        rs[b, 2] ⊻= rs[a, 2] ⊻ (prod_phase_2 != 0)
    end
    return state
end

Base.@propagate_inbounds @inline function measure_z!(state::InvStabilizerState, a;
                                                     force=nothing)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    return _measure_z!(state, n, a, force)
end

# Randomly pick a result
function _measure_z!(state::InvStabilizerState, n, a, force)
    assume(n == state.n)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(n > 0)
    VT2 = Vec{2,ChT}
    lane = VecRange{2}(0)
    local pchunk0, pmask0, zxa
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)
    i0_end = nchunks >> 1
    @inbounds for i0 in 1:i0_end
        i = i0 * 2 - 1
        zxa = xzs[lane + i, 3, a]
        if reduce(|, zxa) == 0
            continue
        end
        pchunk0 = i
        if zxa[1] == 0
            pmask0 = Vec((zero(ChT), one(ChT) << trailing_zeros(zxa[2])))
        else
            pmask0 = Vec((one(ChT) << trailing_zeros(zxa[1]), zero(ChT)))
        end
        @goto rand_measure
    end
    return @inbounds(rs[a, 2]), true

    @label rand_measure
    res = force !== nothing ? force : rand(Bool)

    i0_start = (pchunk0 >> 1) + 2
    if i0_start > i0_end
        @goto single_chunk
    end
    ws = state.ws
    assume(size(ws, 1) == 4)
    chunk_count = 0
    zcum_lo = zero(VT2)
    zcum_hi = zero(VT2)
    @inbounds for i0 in i0_start:i0_end
        i = i0 * 2 - 1
        cnot_tgt = xzs[lane + i, 3, a]
        if reduce(|, cnot_tgt) == 0
            continue
        end
        xzs[lane + i, 3, a] = zero(VT2)

        chunk_count += 1
        ws[lane + 1, chunk_count] = cnot_tgt
        ws[3, chunk_count] = i

        z = xzs[lane + i, 4, a]
        cnot_z = z & cnot_tgt
        zcum_hi ⊻= zcum_lo & cnot_z
        zcum_lo ⊻= cnot_z
    end
    if chunk_count <= 0
        @goto single_chunk
    end
    @inbounds begin
        xzs[lane + pchunk0, 3, a] = zero(VT2)
        cnot_tgt0 = zxa & ~pmask0
        z = xzs[lane + pchunk0, 4, a]
        xzs[lane + pchunk0, 4, a] = z | pmask0
        pz = reduce(|, z & pmask0) != 0
        cnot_z = z & cnot_tgt0
        new_zcum_lo = zcum_lo ⊻ cnot_z
        zcum_lo_cnt_u8 = reduce(+, vcount_ones_u8(new_zcum_lo))
        was_y = pz ⊻ (zcum_lo_cnt_u8 & 1 != 0)

        zcum_hi ⊻= zcum_lo & cnot_z

        zcum_hi_cnt_u8 = reduce(+, vcount_ones_u8(zcum_hi))
        zcum_cnt_u8 = (zcum_hi_cnt_u8 << 1) + zcum_lo_cnt_u8

        cnot_phase = zcum_lo_cnt_u8 ⊻ (ifelse(pz, zcum_cnt_u8,
                                              zcum_cnt_u8 + 0x1) >> 1)
        flip_res = res ⊻ rs[a, 2] ⊻ (cnot_phase & 1 != 0)
        rs[a, 2] = res
    end

    @inbounds for j in 1:n
        for k in 1:2
            if k == 2 && j == a
                continue
            end
            x0 = xzs[lane + pchunk0, 2k - 1, j]
            px = reduce(|, x0 & pmask0) != 0
            z0 = xzs[lane + pchunk0, 2k, j]
            pz = reduce(|, z0 & pmask0) != 0
            cnot_z0 = z0 & cnot_tgt0
            zcum_lo = cnot_z0
            if px
                zcum_hi = zero(VT2)
                xzcum = cnot_z0 & x0
                for cid in 1:chunk_count
                    cnot_tgt = ws[lane + 1, cid]
                    i = ws[3, cid]
                    x = xzs[lane + i, 2k - 1, j]
                    z = xzs[lane + i, 2k, j]
                    cnot_z = z & cnot_tgt
                    xzs[lane + i, 2k - 1, j] = x ⊻ cnot_tgt
                    zcum_hi ⊻= zcum_lo & cnot_z
                    xzcum ⊻= cnot_z & x
                    zcum_lo ⊻= cnot_z
                end
                x0 = x0 ⊻ cnot_tgt0
                xzcum_cnt_u8 = reduce(+, vcount_ones_u8(xzcum))
                zcum_hi_cnt_u8 = reduce(+, vcount_ones_u8(zcum_hi))
                zcum_lo_cnt_u8 = reduce(+, vcount_ones_u8(zcum_lo))
                zcum_cnt_u8 = (zcum_hi_cnt_u8 << 1) + zcum_lo_cnt_u8
                cnot_phase = xzcum_cnt_u8 ⊻ (ifelse(pz, zcum_cnt_u8,
                                                    zcum_cnt_u8 + 0x1) >> 1)
                r = rs[j, k] ⊻ (cnot_phase & 1 != 0)
            else
                for cid in 1:chunk_count
                    cnot_tgt = ws[lane + 1, cid]
                    i = ws[3, cid]
                    z = xzs[lane + i, 2k, j]
                    zcum_lo ⊻= z & cnot_tgt
                end
                zcum_lo_cnt_u8 = reduce(+, vcount_ones_u8(zcum_lo))
                r = rs[j, k]
            end
            new_pz = pz ⊻ (zcum_lo_cnt_u8 & 1 != 0)
            final_pz = ifelse(was_y, px, px ⊻ new_pz)
            final_px = ifelse(was_y, px ⊻ new_pz, new_pz)
            xzs[lane + pchunk0, 2k - 1, j] = _setbit(x0, final_px, pmask0)
            xzs[lane + pchunk0, 2k, j] = _setbit(z0, final_pz, pmask0)
            rs[j, k] = r ⊻ (final_pz & flip_res)
        end
    end
    return res, false

    @label single_chunk

    @inbounds begin
        xzs[lane + pchunk0, 3, a] = zero(VT2)
        cnot_tgt = zxa & ~pmask0
        z = xzs[lane + pchunk0, 4, a]
        xzs[lane + pchunk0, 4, a] = z | pmask0
        pz = reduce(|, z & pmask0) != 0
        if reduce(|, cnot_tgt) != 0
            cnot_z = z & cnot_tgt
            xzcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_z & zxa))
            zcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_z))
            was_y = pz ⊻ (zcum_cnt_u8 & 1 != 0)
            cnot_phase = xzcum_cnt_u8 ⊻ (ifelse(pz, zcum_cnt_u8,
                                                zcum_cnt_u8 + 0x1) >> 1)
            flip_res = res ⊻ rs[a, 2] ⊻ (cnot_phase & 1 != 0)
            rs[a, 2] = res
            for j in 1:n
                for k in 1:2
                    if k == 2 && j == a
                        continue
                    end
                    x = xzs[lane + pchunk0, 2k - 1, j]
                    z = xzs[lane + pchunk0, 2k, j]
                    px = reduce(|, x & pmask0) != 0
                    pz = reduce(|, z & pmask0) != 0

                    cnot_z = z & cnot_tgt
                    zcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_z))
                    new_pz = pz ⊻ (zcum_cnt_u8 & 1 != 0)
                    final_pz = ifelse(was_y, px, px ⊻ new_pz)
                    final_px = ifelse(was_y, px ⊻ new_pz, new_pz)
                    r = rs[j, k] ⊻ (final_pz & flip_res)
                    if px
                        xzcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_z & x))
                        cnot_phase = xzcum_cnt_u8 ⊻ (ifelse(pz, zcum_cnt_u8,
                                                            zcum_cnt_u8 + 0x1) >> 1)
                        x = x ⊻ cnot_tgt
                        r ⊻= cnot_phase & 1 != 0
                    end
                    xzs[lane + pchunk0, 2k - 1, j] = _setbit(x, final_px, pmask0)
                    xzs[lane + pchunk0, 2k, j] = _setbit(z, final_pz, pmask0)
                    rs[j, k] = r
                end
            end
        else
            was_y = pz
            flip_res = res ⊻ rs[a, 2]
            rs[a, 2] = res

            for j in 1:n
                for k in 1:2
                    if k == 2 && j == a
                        continue
                    end
                    x = xzs[lane + pchunk0, 2k - 1, j]
                    z = xzs[lane + pchunk0, 2k, j]
                    px = reduce(|, x & pmask0) != 0
                    pz = reduce(|, z & pmask0) != 0
                    final_pz = ifelse(was_y, px, px ⊻ pz)
                    final_px = ifelse(was_y, px ⊻ pz, pz)
                    xzs[lane + pchunk0, 2k - 1, j] = _setbit(x, final_px, pmask0)
                    xzs[lane + pchunk0, 2k, j] = _setbit(z, final_pz, pmask0)
                    rs[j, k] ⊻= final_pz & flip_res
                end
            end
        end
    end
    return res, false
end

const _StabilizerState = Union{StabilizerState,InvStabilizerState}

@inline init_state_z!(state::_StabilizerState) = init_state_z!(state, false)
@inline init_state_x!(state::_StabilizerState) = init_state_x!(state, false)

function measure_x!(state::_StabilizerState, a; force=nothing)
    @boundscheck check_qubit_bound(state.n, a)
    @inbounds apply!(state, HGate(), a)
    res = @inbounds measure_z!(state, a; force=force)
    @inbounds apply!(state, HGate(), a)
    return res
end

function measure_y!(state::_StabilizerState, a; force=nothing)
    @boundscheck check_qubit_bound(state.n, a)
    @inbounds apply!(state, SXGate(), a)
    res = @inbounds measure_z!(state, a; force=force)
    @inbounds apply!(state, ISXGate(), a)
    return res
end

function measure_zs!(state::_StabilizerState, idxs; force=nothing)
    if isempty(idxs)
        return false, true
    end
    idx0 = @inbounds idxs[1]
    nidxs = length(idxs)
    @boundscheck check_qubit_bound(state.n, idx0)
    if nidxs == 1
        return @inbounds measure_z!(state, idx0; force=force)
    end
    for i in 2:nidxs
        idxi = @inbounds idxs[i]
        @boundscheck check_qubit_bound(state.n, idxi)
        @inbounds apply!(state, CNOTGate(), idxi, idx0)
    end
    res = @inbounds measure_z!(state, idx0; force=force)
    @inbounds for i in 2:nidxs
        apply!(state, CNOTGate(), idxs[i], idx0)
    end
    return res
end

function measure_xs!(state::_StabilizerState, idxs; force=nothing)
    for idx in idxs
        apply!(state, HGate(), idx)
    end
    res = @inbounds measure_zs!(state, idxs; force=force)
    @inbounds for idx in idxs
        apply!(state, HGate(), idx)
    end
    return res
end

function measure_ys!(state::_StabilizerState, idxs; force=nothing)
    for idx in idxs
        apply!(state, SXGate(), idx)
    end
    res = @inbounds measure_zs!(state, idxs; force=force)
    @inbounds for idx in idxs
        apply!(state, ISXGate(), idx)
    end
    return res
end

function measure_paulis!(state::_StabilizerState, xs, zs; force=nothing)
    @assert length(xs) == state.n
    @assert length(zs) == state.n
    local idx0
    @inbounds for i in 1:state.n
        x = xs[i]
        z = zs[i]
        if x
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
            apply!(state, CNOTGate(), i, idx0)
        end
    end
    if !@isdefined(idx0)
        return false, true
    end
    res = @inbounds measure_z!(state, idx0; force=force)
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

function measure_paulis!(state::_StabilizerState, str::PauliString{Bool}; force=nothing)
    if force !== nothing
        force ⊻= str.rs[]
    end
    v, det = measure_paulis!(state, str.xs, str.zs; force=force)
    v ⊻= str.rs[]
    return v, det
end

function measure_stabilizer(state::_StabilizerState, xs, zs, r=false)
    v, det = measure_paulis!(state, xs, zs)
    return v ⊻ r
end
function measure_stabilizer_x(state::_StabilizerState, xs, r=false)
    v, det = measure_xs!(state, xs)
    return v ⊻ r
end
function measure_stabilizer_z(state::_StabilizerState, zs, r=false)
    v, det = measure_zs!(state, zs)
    return v ⊻ r
end

Base.@propagate_inbounds @inline function inject_pauli!(state::_StabilizerState,
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

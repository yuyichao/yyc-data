#!/usr/bin/julia

module Clifford

using Random
using ..Utils
using SIMD

##
# Implement simulation of clifford circuit using the stabilizer state description

##
# Representation of the Pauli string
# The state is generally represented by a collection of Pauli strings.
# Each n-qubit pauli string is described using 2n + 1 bits,
# n X-bits and n Z-bits and one phase bit (R).
# The X and Z bit on each qubit represent the pauli operator on the bit,
#     x=0, z=0: I
#     x=1, z=0: X
#     x=1, z=1: Y
#     x=0, z=1: Z
# and the phase bit R is the sign on the complete Pauli string
# (i.e. a bit 1 means a negative sign).
# For example, an `-XY` pauli string is encoded as `x1=1, z1=0, x2=1, z2=1, r=1`.
# Note that this encoding can only represent Hermitian Pauli strings
# (i.e. without `i` phase) and these are the only ones we need to represent states.

# For most versions of the state representation, we pack these bits into
# larger (UInt64) integers so that we can use bitwise operator to implement
# gate and measurement operations in parallel.

##
# State representations
# There are currently 3 different representations of state.
# 1. Collections of Pauli operators.
#    This can be used to propagate any Pauli operators of interest
#    and is also useful to propagate errors (i.e. diff of the circuits).
# 2. Full set of stabilizer of the state. (CHP)
#    Using the representation out-lined in https://arxiv.org/abs/quant-ph/0406196
# 3. Full set of inverse Pauli frames. (Stim)
#    Following https://arxiv.org/abs/2103.02202
#
# Both the second and the third representations can represent the whole state.
# For most operations including the all operations for large qubit numbers
# the third version (Stim) is faster.
# The second version may be faster for certain operations
# (CNOT and random measurements) on small to intermediate (<1000 qubit) problem sizes.

# Utility functions
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
@inline _getbit(chunk::Vec, mask::Vec) = reduce(|, chunk & mask) != 0
@inline _setbit(chunk, val, mask) = ifelse(val, chunk | mask, chunk & ~mask)

# It turns out that on both x86 and aarch64,
# vectorized popcount is implemented byte-wise
# (the result of the byte-wise popcount is then added together horizontally).
# Since we really only care about the last one or two bits of the result
# casting to `UInt8` vector and do as many operation on it as we can
# reduces the number of horizontal additions we need.
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

# Bool in julia always get masked when loaded even though we know that
# the stored value is either `0` or `1`. Store them using `UInt8` instead
# and do the conversion manually with an explicit assumption of the value range.
@inline function u8_to_bool(v)
    assume(v <= 0x1)
    return v != 0
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

"""
    _combine_1q(x, z, r, xx, xz, xr, zx, zz, zr)

Return the combination of the two Pauli strings X: `xx, xz, xr`
and Z: `zx, zz, zr`, based on the control input `x, z, r`.
"""
@inline function _combine_1q(x, z, r, xx, xz, xr, zx, zz, zr)
    x′ = (x & xx) ⊻ (z & zx)
    z′ = (x & xz) ⊻ (z & zz)
    r ⊻= (z & zr) ⊻ (x & xr) ⊻ ((x & z) & (zz ⊻ xx ⊻ (zx | xz)))
    return x′, z′, r
end
@inline function _anti_commute_1q(x1, z1, x2, z2)
    return (x1 & z2) ⊻ (z1 & x2)
end

# We support the most generic set of single qubit Clifford gates.
"""
    Generic1Q{XX,XZ,XR,ZX,ZZ,ZR}

Generic one qubit Clifford rotation.
The gate is represented by the operation done on the X and the Z operators.
The X operator is turned into the Pauli operator represented by `XX, XZ, XR`
and the Z operator is turned into the Pauli operator represented by `ZX, ZZ, ZR`.

Note that the X and Z Pauli string must anti-commute to maintain the unitarity
of the gate and this is checked in the constructor.
"""
struct Generic1Q{XX,XZ,XR,ZX,ZZ,ZR} <: Clifford1Q
    function Generic1Q{XX,XZ,XR,ZX,ZZ,ZR}() where {XX,XZ,XR,ZX,ZZ,ZR}
        if !_anti_commute_1q(XX, XZ, ZX, ZZ)
            error("X and Z mapping must anti-commute")
        end
        return new{XX,XZ,XR,ZX,ZZ,ZR}()
    end
end
"""
    *(g1::Generic1Q, g2::Generic1Q)

Return the gate corresponding to apply `g1` first and then `g2`.
"""
function Base.:*(::Generic1Q{XX1,XZ1,XR1,ZX1,ZZ1,ZR1},
                 ::Generic1Q{XX2,XZ2,XR2,ZX2,ZZ2,ZR2}) where {XX1,XZ1,XR1,ZX1,ZZ1,ZR1,
                                                              XX2,XZ2,XR2,ZX2,ZZ2,ZR2}
    XX, XZ, XR = _combine_1q(XX1, XZ1, XR1, XX2, XZ2, XR2, ZX2, ZZ2, ZR2)
    ZX, ZZ, ZR = _combine_1q(ZX1, ZZ1, ZR1, XX2, XZ2, XR2, ZX2, ZZ2, ZR2)
    return Generic1Q{XX,XZ,XR,ZX,ZZ,ZR}()
end
"""
    _inv_1q(XX, XZ, XR, ZX, ZZ, ZR)

Return the X and Z Pauli strings representing the inverse of the 1Q gate.
"""
@inline function _inv_1q(XX, XZ, XR, ZX, ZZ, ZR)
    XX′ = ZZ
    XZ′ = XZ
    ZX′ = ZX
    ZZ′ = XX
    XR′ = (ZZ & XR) ⊻ (~(ZR ⊻ ZX) & XZ)
    ZR′ = (XX & ZR) ⊻ (~(XR ⊻ XZ) & ZX)
    return XX′, XZ′, XR′, ZX′, ZZ′, ZR′
end
@inline function Base.inv(::Generic1Q{XX,XZ,XR,ZX,ZZ,ZR}) where {XX,XZ,XR,ZX,ZZ,ZR}
    XX′, XZ′, XR′, ZX′, ZZ′, ZR′ = _inv_1q(XX, XZ, XR, ZX, ZZ, ZR)
    return Generic1Q{XX′,XZ′,XR′,ZX′,ZZ′,ZR′}()
end
const _named_gates = Dict((true, false, false, false, true, false)=>"I",
                          (false, true, false, true, false, false)=>"H",
                          (true, true, false, false, true, true)=>"HXY",
                          (true, false, true, true, true, false)=>"HYZ",
                          (true, false, false, false, true, true)=>"X",
                          (true, false, true, false, true, true)=>"Y",
                          (true, false, true, false, true, false)=>"Z",
                          (true, true, false, false, true, false)=>"S",
                          (true, true, true, false, true, false)=>"IS",
                          (true, false, false, true, true, true)=>"SX",
                          (true, false, false, true, true, false)=>"ISX",
                          (false, true, true, true, false, false)=>"SY",
                          (false, true, false, true, false, true)=>"ISY",
                          (true, true, false, true, false, false)=>"CXYZ",
                          (false, true, false, true, true, false)=>"CZYX")
function _to_pauli_name(x, z, r)
    return (r ? "-" : "+") * (x ? (z ? "Y" : "X") : "Z")
end
function Base.show(io::IO, ::Generic1Q{XX,XZ,XR,ZX,ZZ,ZR}) where {XX,XZ,XR,ZX,ZZ,ZR}
    name = get(_named_gates, (XX,XZ,XR,ZX,ZZ,ZR), nothing)
    if name === nothing
        YX = XX ⊻ ZX
        YZ = XZ ⊻ ZZ
        YR = XR ⊻ ZR ⊻ ZZ ⊻ XX ⊻ (ZX | XZ)
        xname = _to_pauli_name(XX, XZ, XR)
        yname = _to_pauli_name(YX, YZ, YR)
        zname = _to_pauli_name(ZX, ZZ, ZR)
        print(io, "Generic1Q(I->+I, X->$(xname), Y->$(yname), Z->$(zname))")
    else
        print(io, "$(name)Gate()")
    end
end
# Create convinience alias for certain known-named gates.
for (param, name) in _named_gates
    Generic1Q{param...}() # Builtin test
    @eval const $(Symbol("$(name)Gate")) = $(Generic1Q{param...})
end
const SZGate = SGate
const ISZGate = ISGate
const HXZGate = HGate
const HZXGate = HGate
const HYXGate = HXYGate
const HZYGate = HYZGate

# Implementation of forward-propagation of Pauli
Base.@propagate_inbounds @inline function apply!(::Generic1Q{XX,XZ,XR,ZX,ZZ,ZR},
                                                 xas, zas, rs) where {XX,XZ,XR,ZX,ZZ,ZR}
    YPHASE = ZZ ⊻ XX ⊻ (ZX | XZ)
    NEED_X = ZX | XR | YPHASE
    NEED_Z = XZ | ZR | YPHASE
    NEED_R = XR | ZR | YPHASE

    xa = NEED_X ? xas[] : zero(eltype(xas))
    za = NEED_Z ? zas[] : zero(eltype(zas))
    if ZX
        if XX
            xas[] = xa ⊻ za
        else
            xas[] = za
        end
    end
    if XZ
        if ZZ
            zas[] = xa ⊻ za
        else
            zas[] = xa
        end
    end

    if NEED_R
        r = rs[]
        if XR
            r ⊻= xa
        end
        if ZR
            r ⊻= za
        end
        if YPHASE
            r ⊻= xa & za
        end
        rs[] = r
    end
    return
end

## Two qubit gates
abstract type Clifford2Q end
@inline function _anti_commute_2q(x11, z11, x21, z21,
                                  x12, z12, x22, z22)
    ac1 = _anti_commute_1q(x11, z11, x12, z12)
    ac2 = _anti_commute_1q(x21, z21, x22, z22)
    return ac1 ⊻ ac2
end
@inline function _phase_2q_anti_commute(x11, z11, x21, z21,
                                        x12, z12, x22, z22)
    if _anti_commute_1q(x11, z11, x12, z12)
        return z12 ⊻ x11 ⊻ (x12 | z11)
    else
        return z22 ⊻ x21 ⊻ (x22 | z21)
    end
end
@inline function _phase_2q_commute(x11, z11, x21, z21,
                                   x12, z12, x22, z22)
    if _anti_commute_1q(x11, z11, x12, z12)
        return (z12 ⊻ x11 ⊻ (x12 | z11)) & (z22 ⊻ x21 ⊻ (x22 | z21))
    else
        return zero(x11)
    end
end
# For anti commuting string, the phase is computed relative to the phase of (X * Z)
# For commuting string, the phase is absolute.
@inline function _phase_2q(x11, z11, x21, z21,
                           x12, z12, x22, z22)
    if _anti_commute_2q(x11, z11, x21, z21, x12, z12, x22, z22)
        return _phase_2q_anti_commute(x11, z11, x21, z21, x12, z12, x22, z22)
    else
        return _phase_2q_commute(x11, z11, x21, z21, x12, z12, x22, z22)
    end
end
@inline function _combine_2q_nophase(b1, b2,
                                     x11, z11, x21, z21, r1,
                                     x12, z12, x22, z22, r2)
    x1′ = (b1 & x11) ⊻ (b2 & x12)
    z1′ = (b1 & z11) ⊻ (b2 & z12)
    x2′ = (b1 & x21) ⊻ (b2 & x22)
    z2′ = (b1 & z21) ⊻ (b2 & z22)
    r = (b2 & r2) ⊻ (b1 & r1)
    return x1′, z1′, x2′, z2′, r
end
@inline function _combine4_2q(b1, b2, b3, b4, r,
                              x11, z11, x21, z21, r1,
                              x12, z12, x22, z22, r2,
                              x13, z13, x23, z23, r3,
                              x14, z14, x24, z24, r4)

    x1_12, z1_12, x2_12, z2_12, r_12 = _combine_2q_nophase(b1, b2,
                                                           x11, z11, x21, z21, r1,
                                                           x12, z12, x22, z22, r2)
    r_12 ⊻= b1 & b2 & _phase_2q(x11, z11, x21, z21, x12, z12, x22, z22)
    x1_34, z1_34, x2_34, z2_34, r_34 = _combine_2q_nophase(b3, b4,
                                                           x13, z13, x23, z23, r3,
                                                           x14, z14, x24, z24, r4)
    r_34 ⊻= b3 & b4 & _phase_2q(x13, z13, x23, z23, x14, z14, x24, z24)

    x1′ = x1_12 ⊻ x1_34
    z1′ = z1_12 ⊻ z1_34
    x2′ = x2_12 ⊻ x2_34
    z2′ = z2_12 ⊻ z2_34
    r ⊻= r_12 ⊻ r_34
    r ⊻= _phase_2q(x1_12, z1_12, x2_12, z2_12, x1_34, z1_34, x2_34, z2_34)

    return x1′, z1′, x2′, z2′
end

struct CNOTGate <: Clifford2Q
end
# Implementation of forward-propagation of Pauli
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
    assume(length(str.xs) == length(str.zs))
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

# Injecting an Pauli operator
# This is intended for injecting or correcting Pauli errors
# into the state diff.
@inline function inject_pauli!(str::PauliString, x, z, i)
    @boundscheck check_qubit_bound(str.n, i)
    @inbounds str.xs[i] ⊻= x
    @inbounds str.zs[i] ⊻= z
    # Ignore signs
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

# Implementation of stablizer state based on https://arxiv.org/abs/quant-ph/0406196
# The bit in each Pauli strings corresponds to each qubit is packed into a `UInt64`
# so that we can do parallel gate operations on multiple strings at the same time.
# We store in total `2n` strings `n` for the destablizers and stabilizers each.
# The X and Z terms for each string are stored in two matrices separately.
struct StabilizerState
    n::Int
    xs::Matrix{ChT}
    zs::Matrix{ChT}
    rs::Matrix{ChT}
    # This is the 2n+1 string in the CHP paper as workspace
    # to figure out the result of deterministic measurements.
    # This relaxes some requirement during measurement allowing us to
    # do out-of-bound read/write on the unused bits at the end
    # without having to mask them out.
    wxzs::Matrix{ChT}
    function StabilizerState(n)
        nchunks = (2n - 1) ÷ _chunk_len + 1
        xs = zeros(ChT, nchunks, n)
        zs = zeros(ChT, nchunks, n)
        # Use a matrix here to help the julia compiler realizing
        # that this cannot be resized and that the size and pointer
        # are both lifetime constants.
        rs = zeros(ChT, nchunks, 1)
        assume(size(xs, 1) == size(zs, 1) == size(rs, 1))
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
    nchunks = size(rs, 1)
    assume(nchunks > 0)
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1))
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
    nchunks = size(rs, 1)
    assume(nchunks > 0)
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1))
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

# Initialize all qubits to the eigenstates of `Z` with sign `v`
# (true for negative sign)
function init_state_z!(state::StabilizerState, v::Bool)
    n = state.n
    xs = state.xs
    zs = state.zs
    rs = state.rs
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1) == length(rs))
    assume(length(xs) == length(zs))
    xs .= 0
    zs .= 0
    assume(n > 0)
    @inbounds for i in 1:n
        # Destablizer to X
        chunk1, mask1 = _get_chunk_mask(i)
        xs[chunk1, i] = mask1
        # Stablizer to Z
        chunk2, mask2 = _get_chunk_mask(i + n)
        zs[chunk2, i] = mask2
    end
    rs .= _cast_bits(ChT, v)
    return state
end

# Initialize all qubits to the eigenstates of `X` with sign `v`
# (true for negative sign)
function init_state_x!(state::StabilizerState, v::Bool)
    n = state.n
    xs = state.xs
    zs = state.zs
    rs = state.rs
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1) == length(rs))
    assume(length(xs) == length(zs))
    xs .= 0
    zs .= 0
    assume(n > 0)
    @inbounds for i in 1:n
        # Destablizer to Z
        chunk1, mask1 = _get_chunk_mask(i)
        zs[chunk1, i] = mask1
        # Stablizer to X
        chunk2, mask2 = _get_chunk_mask(i + n)
        xs[chunk2, i] = mask2
    end
    rs .= _cast_bits(ChT, v)
    return state
end

# Initialize all qubits to the eigenstates of `Y` with sign `v`
# (true for negative sign)
function init_state_y!(state::StabilizerState, v::Bool)
    n = state.n
    xs = state.xs
    zs = state.zs
    rs = state.rs
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1) == length(rs))
    assume(length(xs) == length(zs))
    xs .= 0
    zs .= 0
    assume(n > 0)
    @inbounds for i in 1:n
        # Destablizer to X
        chunk1, mask1 = _get_chunk_mask(i)
        # Stablizer to Y
        chunk2, mask2 = _get_chunk_mask(i + n)
        if chunk1 == chunk2
            xs[chunk1, i] = mask1 | mask2
        else
            xs[chunk1, i] = mask1
            xs[chunk2, i] = mask2
        end
        zs[chunk2, i] = mask2
    end
    rs .= _cast_bits(ChT, v)
    return state
end

# Accumulate the phase (±i) of the product of the two Pauli operators
# `x1, z1` and `x2, z2` into the two-bit counter `hi, lo`.
@inline function _accum_pauli_prod_phase(hi, lo, x1, z1, x2, z2)
    v1 = x1 & z2
    v2 = x2 & z1
    # This is different from the efficient implementation used in Stim.
    # The Stim expression (used below) is more efficient when the
    # new value of x and z (i.e. `x1 ⊻ x2` and `z1 ⊻ z2`) are already computed.
    # The following expression is more efficient without those pre-computed values
    # and is optimized by enumeration.
    m = z1 ⊻ x2 ⊻ (x1 | z2)
    change = v1 ⊻ v2
    hi = hi ⊻ ((m ⊻ lo) & change)
    lo = lo ⊻ change
    return hi, lo
end

# These are implementations of the `rowsum` function from CHP,
# which is basically a multiplication between communing Pauli strings.
# The different versions here corresponds to different usage patterns
# and different ways for parallel processing.
# All of these functions have a single user in the `_measure_z!`
# for `StabilizerState` and should all be inlined into that function.

# Multiply the single Pauli string located at `chunk_i, mask_i`
# onto the Pauli strings at `chunk_h, mask_h`.
# (i.e. `mask_i` have only a single bit set whereas `mask_h` might have multiple set)
@inline function _pauli_rowsum1!(state::StabilizerState, chunk_h, mask_h,
                                 chunk_i, mask_i)
    xs = state.xs
    zs = state.zs
    rs = state.rs
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1))
    hi = zero(ChT)
    lo = zero(ChT)
    assume(state.n > 0)
    @inbounds for j in 1:state.n
        xi = _getbit(xs[chunk_i, j], mask_i)
        zi = _getbit(zs[chunk_i, j], mask_i)
        # Based on testing, it seems that the cost of the branch is less than
        # the saving in other operations (probably memory read/writes)
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
        rh = rs[chunk_h]
        hi = ifelse(_getbit(rs[chunk_i], mask_i), ~hi, hi)
        rs[chunk_h] = rh ⊻ (hi & mask_h)
    end
    return state
end

# Similar to `_pauli_rowsum1!` but without branches.
# This seems to perform better for small qubit count.
@inline function _pauli_rowsum1_2!(state::StabilizerState, chunk_h, mask_h,
                                   chunk_i, mask_i)
    xs = state.xs
    zs = state.zs
    rs = state.rs
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1))
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
        rh = rs[chunk_h]
        hi = ifelse(_getbit(rs[chunk_i], mask_i), ~hi, hi)
        rs[chunk_h] = rh ⊻ (hi & mask_h)
    end
    return state
end

# Multiply pair-wise the Pauli strings in `chunk_i, mask_i` onto the Pauli strings
# in the workspace. The workspace phase is kept track of in `wrs`.
@inline function _pauli_rowsum2!(state::StabilizerState, chunk_i, mask_i, wrs)
    n = state.n
    xs = state.xs
    zs = state.zs
    rs = state.rs
    wxzs = state.wxzs
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1))
    assume(size(wxzs, 1) == n)
    hi = zero(ChT)
    lo = zero(ChT)
    assume(n > 0)
    @inbounds for j in 1:n
        xi = xs[chunk_i, j]
        xh = wxzs[j, 1]
        zi = zs[chunk_i, j]
        zh = wxzs[j, 2]
        hi, lo = _accum_pauli_prod_phase(hi, lo, xi, zi, xh, zh)
        wxzs[j, 1] = xh ⊻ (xi & mask_i)
        wxzs[j, 2] = zh ⊻ (zi & mask_i)
    end
    @inbounds begin
        wrs[] ⊻= (hi ⊻ rs[chunk_i, 1]) & mask_i
    end
    return state
end

# When there are multiple non-trivial Pauli string in the workspace,
# Compute the phase when all of them are multiplied together.
# We do not care about the string itself in this case.
@inline function _pauli_rowsum3!(state::StabilizerState, wrs)
    n = state.n
    @assert sizeof(ChT) == 8
    wxzs = state.wxzs
    hi = zero(ChT)
    lo = zero(ChT)
    assume(size(wxzs, 1) == n)
    assume(n > 0)
    @inbounds for j in 1:n
        # Although we can assume all of the Pauli strings commutes with each other
        # therefore the phase should be independent of the multiplication order,
        # I cannot find an algorithm that is explicit order-independent.
        # As such, the way I'm computing this is to simulate a log(64) steps
        # folding reduction. However, instead of computing the phase
        # for each multiplication I'm computing and recording the operand
        # of each mulplication (there are in total 63 of them) in two vectors
        # (`xh, zh` and `xi, zi`) and then computing the sum of the phase accumulated
        # on all of the 63 multiplications.

        x = wxzs[j, 1]
        z = wxzs[j, 2]

        # Multiply 64 bits into 32 bits
        x1 = (x >> 32) % UInt32
        z1 = (z >> 32) % UInt32
        x2 = x % UInt32
        z2 = z % UInt32

        # Record the operands in the lowest 32 bits [31...0]
        xh = x1 % UInt64
        zh = z1 % UInt64
        xi = x2 % UInt64
        zi = z2 % UInt64

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        # Multiply 32 bits into 16 bits
        x1 = x >> 16
        z1 = z >> 16
        x2 = x & 0x0000ffff
        z2 = z & 0x0000ffff

        # Record the operands in bits [47...32]
        xh |= (x1 % UInt64) << 32
        zh |= (z1 % UInt64) << 32
        xi |= (x2 % UInt64) << 32
        zi |= (z2 % UInt64) << 32

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        # Multiply 16 bits into 8 bits
        x1 = x >> 8
        z1 = z >> 8
        x2 = x & 0x000000ff
        z2 = z & 0x000000ff

        # Record the operands in bits [55...48]
        xh |= (x1 % UInt64) << 48
        zh |= (z1 % UInt64) << 48
        xi |= (x2 % UInt64) << 48
        zi |= (z2 % UInt64) << 48

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        # Multiply 8 bits into 4 bits
        x1 = x >> 4
        z1 = z >> 4
        x2 = x & 0x0000000f
        z2 = z & 0x0000000f

        # Record the operands in bits [59...56]
        xh |= (x1 % UInt64) << 56
        zh |= (z1 % UInt64) << 56
        xi |= (x2 % UInt64) << 56
        zi |= (z2 % UInt64) << 56

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        # Multiply 4 bits into 2 bits
        x1 = x >> 2
        z1 = z >> 2
        x2 = x & 0x00000003
        z2 = z & 0x00000003

        # Record the operands in bits [61...60]
        xh |= (x1 % UInt64) << 60
        zh |= (z1 % UInt64) << 60
        xi |= (x2 % UInt64) << 60
        zi |= (z2 % UInt64) << 60

        x = x1 ⊻ x2
        z = z1 ⊻ z2

        # Multiply 2 bits into 1 bits
        x1 = x >> 1
        z1 = z >> 1
        x2 = x & 0x00000001
        z2 = z & 0x00000001

        # Record the operands in bit 62
        xh |= (x1 % UInt64) << 62
        zh |= (z1 % UInt64) << 62
        xi |= (x2 % UInt64) << 62
        zi |= (z2 % UInt64) << 62

        hi, lo = _accum_pauli_prod_phase(hi, lo, xh, zh, xi, zi)
    end
    return (count_ones(wrs[] ⊻ hi) ⊻ (count_ones(lo) >> 1)) & 1 != 0
end

# With the rearrangement of the stabilizer finished, set the (p-n)-th stabilizer
# to reflect the measurement result. `chunk2, mask2` is the bit index/mask
# for the stabilizer (i.e. the p-th string)
@inline function _inject_zresult!(state::StabilizerState, n, a, p, chunk2, mask2, res)
    chunk1, mask1 = _get_chunk_mask(p - n) # index and mask for the de-stabilizer
    xs = state.xs
    zs = state.zs
    rs = state.rs
    assume(n > 0)
    assume(size(xs, 1) == size(zs, 1) == size(rs, 1))
    @inbounds for i in 1:n
        # Set the destabilizer to the current stabilizer
        # and set the stabilizer to the correct `Z` with the correct sign.

        x2 = xs[chunk2, i]
        xs[chunk2, i] = x2 & ~mask2
        # Note that the load of x1 has to happen after the store of x2
        # since the two may alias
        xs[chunk1, i] = _setbit(xs[chunk1, i], _getbit(x2, mask2), mask1)

        z2 = zs[chunk2, i]
        zs[chunk2, i] = _setbit(z2, i == a, mask2)
        # Note that the load of z1 has to happen after the store of z2
        # since the two may alias
        zs[chunk1, i] = _setbit(zs[chunk1, i], _getbit(z2, mask2), mask1)
    end
    @inbounds begin
        r2 = rs[chunk2, 1]
        rs[chunk2, 1] = _setbit(r2, res, mask2)
        # Note that the load of r1 has to happen after the store of r2
        # since the two may alias
        rs[chunk1, 1] = _setbit(rs[chunk1, 1], _getbit(r2, mask2), mask1)
    end
    return
end

Base.@propagate_inbounds @inline function measure_z!(state::StabilizerState, a;
                                                     force=nothing)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    return _measure_z!(state, n, a, force)
end

function _measure_z!(state::StabilizerState, n, a, force)
    assume(n == state.n)
    p = 0
    chunk_p = 0
    mask_p = zero(ChT)
    found_p = false
    assume(n > 0)
    xs = state.xs
    rs = state.rs
    nchunks = size(rs, 1)
    assume(nchunks > 0)
    assume(size(xs, 1) == size(rs, 1))
    @inbounds for _p in (n + 1):(2 * n)
        chunk_p, mask_p = _get_chunk_mask(_p)
        if _getbit(xs[chunk_p, a], mask_p)
            p = _p
            found_p = true
            break
        end
    end
    @inbounds if found_p
        # Rearrange the stabilizers so that there's a single stablizer that is not
        # communting with our measurement.
        if nchunks <= 1
            # Special case optimizing for small bit count
            mask_i = xs[1, a]
            bitidx = (p - 1) % UInt
            mask_i = ifelse(bitidx < _chunk_len,
                            mask_i & ~(one(ChT) << bitidx), mask_i)
            if mask_i != 0
                _pauli_rowsum1_2!(state, 1, mask_i, chunk_p, mask_p)
            end
        else
            for i in 1:nchunks
                mask_i = xs[i, a]
                bitidx = (p - 1 - (i - 1) * _chunk_len) % UInt
                mask_i = ifelse(bitidx < _chunk_len,
                                mask_i & ~(one(ChT) << bitidx), mask_i)
                if mask_i == 0
                    continue
                end
                _pauli_rowsum1!(state, i, mask_i, chunk_p, mask_p)
            end
        end
        # Randomly pick a result
        res = force !== nothing ? force : rand(Bool)
        _inject_zresult!(state, n, a, p, chunk_p, mask_p, res)
        return res, false
    else
        # Multiply all the stabilizers that are needed to construct our
        # measured Z, using the coefficients from the de-stabilizers.
        nfullchunks, rembits = divrem(n, _chunk_len)
        wxzs = state.wxzs
        assume(size(wxzs, 1) == n)
        assume(size(wxzs, 2) == 2)
        assume(length(wxzs) == 2 * n)
        wxzs .= zero(ChT)
        wrs = Ref(zero(ChT))
        total_mask = zero(ChT)
        if rembits == 0
            # If the bit number aligns with the chunk size,
            # then the bit mask for stablizer and destabilizer aligns
            # ans we can simply use the corresponding integers with matching index.
            assume(nfullchunks > 0)
            for i in 1:nfullchunks
                mask_i = xs[i, a]
                if mask_i != 0
                    total_mask |= mask_i
                    _pauli_rowsum2!(state, i + nfullchunks, mask_i, wrs)
                end
            end
        else
            # If n is not aligned to a full chunk, we need to shift the destablizer
            # mask to match the stablizer mask.
            # As we iterate, keep track of the destabilizer mask that we've loade
            # but have not used yet to be used on the next iteration.
            mask = zero(ChT)
            for i in 1:nfullchunks
                new_mask = xs[i, a]
                mask_i = mask | (new_mask << rembits)
                mask = new_mask >> (_chunk_len - rembits)
                if mask_i != 0
                    total_mask |= mask_i
                    _pauli_rowsum2!(state, i + nfullchunks, mask_i, wrs)
                end
            end
            new_mask = xs[nfullchunks + 1, a]
            mask_i = mask | (new_mask << rembits)
            rembits2 = rembits * 2
            mask_i &= (one(ChT) << rembits2) - one(ChT)
            if mask_i != 0
                total_mask |= mask_i
                _pauli_rowsum2!(state, 2 * nfullchunks + 1, mask_i, wrs)
            end
            if rembits2 > _chunk_len
                mask_i = new_mask >> (_chunk_len - rembits)
                # Mask out the end of the range
                mask_i &= (one(ChT) << (rembits2 - _chunk_len)) - one(ChT)
                if mask_i != 0
                    total_mask |= mask_i
                    _pauli_rowsum2!(state, 2 * nfullchunks + 2, mask_i, wrs)
                end
            end
        end
        # Mostly to catch the cases where there actually is only a single term
        # but could also catch when we are lucky and the multiple terms aligns.
        if count_ones(total_mask) == 1
            return wrs[] != 0, true
        end
        return _pauli_rowsum3!(state, wrs), true
    end
end

# Implement the inverted stabilizer tableau based on https://arxiv.org/abs/2103.02202.
struct InvStabilizerState
    n::Int
    # The last dimension of `xzs` is the latest time bit index.
    # Each matrix in the first two dimensions stores the beginning-of-time
    # Pauli string that corresponds to `X` and `Z` on this bit at latest time.
    # The signes for the `X` and `Z` Pauli strings are stored
    # in `rs` as an `n` by `2` matrix.
    xzs::Array{ChT,3}
    rs::Matrix{UInt8}
    # This (`ws`) is a workspace member currently used during random measurement
    # to record to beginning-of-time CNOT that needs to be applied.
    ws::Matrix{ChT}
    function InvStabilizerState(n)
        # Align the chunk number to 2 for better alignment and easier processing.
        nchunks = ((n - 1) ÷ (_chunk_len * 2) + 1) * 2
        # XX,            XZ,            ZX,            ZZ
        # X stab X term, X stab Z term, Z stab X term, Z stab Z term
        xzs = zeros(ChT, nchunks, 4, n)
        rs = zeros(UInt8, n, 2)
        assume(size(xzs, 2) == 4)
        @inbounds for i in 1:n
            chunk, mask = _get_chunk_mask(i)
            xzs[chunk, 1, i] = mask
            xzs[chunk, 4, i] = mask
        end
        return new(n, xzs, rs, zeros(ChT, 4, nchunks >> 1))
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
                       u8_to_bool(state.rs[i, k]))
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

function init_state_y!(state::InvStabilizerState, v::Bool)
    n = state.n
    xzs = state.xzs
    assume(size(xzs, 2) == 4)
    xzs .= 0
    assume(n > 0)
    @inbounds for i in 1:n
        chunk, mask = _get_chunk_mask(i)
        xzs[chunk, 1, i] = mask
        xzs[chunk, 3, i] = mask
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

@inline function _pauli_multiply_kernel_8!(px1s, pz1s, px2s, pz2s, i, chi, clo,
                                           ::Val{do_swap}) where {do_swap}
    VT8 = Vec{8,ChT}

    x1 = vload(VT8, _gep_array(px1s, (), (i,)))
    x2 = vload(VT8, _gep_array(px2s, (), (i,)))
    new_x1 = x1 ⊻ x2
    vstore(new_x1, _gep_array(px1s, (), (i,)))
    if do_swap
        vstore(x1, _gep_array(px2s, (), (i,)))
    end
    z1 = vload(VT8, _gep_array(pz1s, (), (i,)))
    z2 = vload(VT8, _gep_array(pz2s, (), (i,)))
    new_z1 = z1 ⊻ z2
    vstore(new_z1, _gep_array(pz1s, (), (i,)))
    if do_swap
        vstore(z1, _gep_array(pz2s, (), (i,)))
    end

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2

    m_1 = Vec((m[1], m[2]))
    m_2 = Vec((m[3], m[4]))
    m_3 = Vec((m[5], m[6]))
    m_4 = Vec((m[7], m[8]))
    change_1 = Vec((change[1], change[2]))
    change_2 = Vec((change[3], change[4]))
    change_3 = Vec((change[5], change[6]))
    change_4 = Vec((change[7], change[8]))

    chi ⊻= (m_1 ⊻ clo) & change_1
    clo ⊻= change_1
    chi ⊻= (m_2 ⊻ clo) & change_2
    clo ⊻= change_2
    chi ⊻= (m_3 ⊻ clo) & change_3
    clo ⊻= change_3
    chi ⊻= (m_4 ⊻ clo) & change_4
    clo ⊻= change_4
    return chi, clo
end

@inline function _pauli_multiply_kernel_4!(px1s, pz1s, px2s, pz2s, chi, clo,
                                           ::Val{do_swap}) where {do_swap}
    VT4 = Vec{4,ChT}

    x1 = vload(VT4, px1s)
    x2 = vload(VT4, px2s)
    new_x1 = x1 ⊻ x2
    vstore(new_x1, px1s)
    if do_swap
        vstore(x1, px2s)
    end
    z1 = vload(VT4, pz1s)
    z2 = vload(VT4, pz2s)
    new_z1 = z1 ⊻ z2
    vstore(new_z1, pz1s)
    if do_swap
        vstore(z1, pz2s)
    end

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    m_1 = Vec((m[1], m[2]))
    m_2 = Vec((m[3], m[4]))
    change_1 = Vec((change[1], change[2]))
    change_2 = Vec((change[3], change[4]))
    if chi !== nothing
        chi ⊻= (m_1 ⊻ clo) & change_1
        clo ⊻= change_1
        chi ⊻= (m_2 ⊻ clo) & change_2
        clo ⊻= change_2
    else
        chi = m_1 & change_1
        clo = change_1
        chi ⊻= (m_2 ⊻ clo) & change_2
        clo ⊻= change_2
    end
    px1s += sizeof(VT4)
    pz1s += sizeof(VT4)
    px2s += sizeof(VT4)
    pz2s += sizeof(VT4)
    return px1s, pz1s, px2s, pz2s, chi, clo
end

@inline function _pauli_multiply_kernel_2!(px1s, pz1s, px2s, pz2s, chi, clo,
                                           ::Val{do_swap}) where {do_swap}
    VT2 = Vec{2,ChT}
    x1 = vloada(VT2, px1s)
    x2 = vloada(VT2, px2s)
    new_x1 = x1 ⊻ x2
    vstorea(new_x1, px1s)
    if do_swap
        vstorea(x1, px2s)
    end
    z1 = vloada(VT2, pz1s)
    z2 = vloada(VT2, pz2s)
    new_z1 = z1 ⊻ z2
    vstorea(new_z1, pz1s)
    if do_swap
        vstorea(z1, pz2s)
    end

    v1 = x1 & z2
    v2 = x2 & z1
    m = new_x1 ⊻ new_z1 ⊻ v1
    change = v1 ⊻ v2
    return chi ⊻ ((m ⊻ clo) & change), clo ⊻ change
end

@inline function _pauli_multiply_kernel_single_2!(px1s, pz1s, px2s, pz2s,
                                                  ::Val{do_swap}) where {do_swap}
    VT4 = Vec{4,ChT}
    xz1 = vload(VT4, px1s)
    xz2 = vload(VT4, px2s)
    new_xz1 = xz1 ⊻ xz2
    vstore(new_xz1, px1s)
    if do_swap
        vstorea(xz1, px2s)
    end

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
    return m & change, change
end

@inline function _pauli_multiply_phase(chi, clo)
    cnt = vcount_ones_u8(clo) + vcount_ones_u8(chi) << 1
    return reduce(+, cnt) & 0x3
end

# P1 = P1 * P2
@inline function pauli_multiply!(px1s, pz1s, px2s, pz2s, n,
                                 ds::Val{do_swap}=Val(false)) where {do_swap}
    VT2 = Vec{2,ChT}

    # Manual jump threading for the small n cases
    if n <= 2
        chi, clo = _pauli_multiply_kernel_single_2!(px1s, pz1s, px2s, pz2s, ds)
        @goto n2_case_end
    elseif n < 8
        px1s, pz1s, px2s, pz2s, chi, clo =
            _pauli_multiply_kernel_4!(px1s, pz1s, px2s, pz2s, nothing, nothing, ds)
        @goto n4_case_end
    end
    chi = zero(VT2)
    clo = zero(VT2)
    nalign = n & ~7
    # LLVM may think the vstore in the loop aliases the array pointer
    # so we extract the array pointer out of the loop.
    @inbounds for i0 in 1:(n >> 3)
        i = (i0 - 1) * 8 + 1
        chi, clo = _pauli_multiply_kernel_8!(px1s, pz1s, px2s, pz2s, i, chi, clo, ds)
    end
    px1s += nalign * 8
    pz1s += nalign * 8
    px2s += nalign * 8
    pz2s += nalign * 8

    if n & 4 == 0
        @goto n4_case_end
    end
    @label n4_case
    px1s, pz1s, px2s, pz2s, chi, clo =
        _pauli_multiply_kernel_4!(px1s, pz1s, px2s, pz2s, chi, clo, ds)
    @label n4_case_end

    if n & 2 == 0
        @goto n2_case_end
    end
    @label n2_case
    chi, clo = _pauli_multiply_kernel_2!(px1s, pz1s, px2s, pz2s, chi, clo, ds)
    @label n2_case_end

    return _pauli_multiply_phase(chi, clo)
end

@inline function pauli_multiply_2!(px1s_1, pz1s_1, px2s_1, pz2s_1,
                                   px1s_2, pz1s_2, px2s_2, pz2s_2, n)
    VT2 = Vec{2,ChT}
    ds = Val(false)

    # Manual jump threading for the small n cases
    if n <= 2
        chi_1, clo_1 = _pauli_multiply_kernel_single_2!(px1s_1, pz1s_1,
                                                        px2s_1, pz2s_1, ds)
        chi_2, clo_2 = _pauli_multiply_kernel_single_2!(px1s_2, pz1s_2,
                                                        px2s_2, pz2s_2, ds)
        @goto n2_case_end
    elseif n < 8
        px1s_1, pz1s_1, px2s_1, pz2s_1, chi_1, clo_1 =
            _pauli_multiply_kernel_4!(px1s_1, pz1s_1, px2s_1, pz2s_1,
                                      nothing, nothing, ds)
        px1s_2, pz1s_2, px2s_2, pz2s_2, chi_2, clo_2 =
            _pauli_multiply_kernel_4!(px1s_2, pz1s_2, px2s_2, pz2s_2,
                                      nothing, nothing, ds)
        @goto n4_case_end
    end
    chi_1 = zero(VT2)
    clo_1 = zero(VT2)
    chi_2 = zero(VT2)
    clo_2 = zero(VT2)
    nalign = n & ~7
    # LLVM may think the vstore in the loop aliases the array pointer
    # so we extract the array pointer out of the loop.
    @inbounds for i0 in 1:(n >> 3)
        i = (i0 - 1) * 8 + 1
        chi_1, clo_1 = _pauli_multiply_kernel_8!(px1s_1, pz1s_1, px2s_1, pz2s_1,
                                                 i, chi_1, clo_1, ds)
        chi_2, clo_2 = _pauli_multiply_kernel_8!(px1s_2, pz1s_2, px2s_2, pz2s_2,
                                                 i, chi_2, clo_2, ds)
    end

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
    px1s_1, pz1s_1, px2s_1, pz2s_1, chi_1, clo_1 =
        _pauli_multiply_kernel_4!(px1s_1, pz1s_1, px2s_1, pz2s_1, chi_1, clo_1, ds)
    px1s_2, pz1s_2, px2s_2, pz2s_2, chi_2, clo_2 =
        _pauli_multiply_kernel_4!(px1s_2, pz1s_2, px2s_2, pz2s_2, chi_2, clo_2, ds)
    @label n4_case_end

    if n & 2 == 0
        @goto n2_case_end
    end
    @label n2_case
    chi_1, clo_1 = _pauli_multiply_kernel_2!(px1s_1, pz1s_1, px2s_1, pz2s_1,
                                             chi_1, clo_1, ds)
    chi_2, clo_2 = _pauli_multiply_kernel_2!(px1s_2, pz1s_2, px2s_2, pz2s_2,
                                             chi_2, clo_2, ds)
    @label n2_case_end

    return _pauli_multiply_phase(chi_1, clo_1), _pauli_multiply_phase(chi_2, clo_2)
end

Base.@propagate_inbounds @inline function apply!(
    state::InvStabilizerState, ::Generic1Q{XX,XZ,XR,ZX,ZZ,ZR},
    a) where {XX,XZ,XR,ZX,ZZ,ZR}

    n = state.n
    @boundscheck check_qubit_bound(n, a)
    xzs = state.xzs
    rs = state.rs
    nchunks = size(xzs, 1)
    assume(nchunks & 1 == 0)
    assume(nchunks >= 2)
    assume(size(xzs, 2) == 4)
    assume(size(rs, 1) == n)

    # 3 different cases
    # 1. no XZ change, only phase change: we need just need to apply sign changes
    # 2. XZ exchange, swap X and Z and apply sign changes
    # 3. one and only one of X and Z was mapped to the original Y
    #    we need to multiply the X and Z strings together
    #    3.1. The other one remains unchanged (up to a sign) (e.g. X <- X, Z <- Y)
    #    3.2. The other one is swapped (e.g. X <- Z, Z <- Y)

    @inbounds if XX & !XZ & !ZX & ZZ
        # Case 1: no XZ change
        if XR
            rs[a, 1] = ~u8_to_bool(rs[a, 1])
        end
        if ZR
            rs[a, 2] = ~u8_to_bool(rs[a, 2])
        end
    elseif !XX & XZ & ZX & !ZZ
        # Case 2: XZ exchange
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
        rx = u8_to_bool(rs[a, 1])
        rz = u8_to_bool(rs[a, 2])
        rs[a, 2] = XR ? ~rx : rx
        rs[a, 1] = ZR ? ~rz : rz
    else
        inv_XX, inv_XZ, inv_XR, inv_ZX, inv_ZZ, inv_ZR =
            _inv_1q(XX, XZ, XR, ZX, ZZ, ZR)
        x_was_y = inv_XX & inv_XZ
        # Both of the cross terms are true, meaning none of X or Z maps to themselves
        do_swap = XZ & ZX

        GC.@preserve xzs begin
            px1s = pointer(@view(xzs[1, x_was_y ? 1 : 3, a]))
            pz1s = pointer(@view(xzs[1, x_was_y ? 2 : 4, a]))
            px2s = pointer(@view(xzs[1, x_was_y ? 3 : 1, a]))
            pz2s = pointer(@view(xzs[1, x_was_y ? 4 : 2, a]))
            prod_phase = pauli_multiply!(px1s, pz1s, px2s, pz2s, nchunks, Val(do_swap))
            assume(prod_phase & 0x1 != 0)
            pr = prod_phase & 0x2 != 0
            rx = u8_to_bool(rs[a, 1])
            rz = u8_to_bool(rs[a, 2])
            if x_was_y
                # Z is multiplied into X, so is the phase
                new_rx = rx ⊻ rz ⊻ pr
                if inv_XR
                    new_rx = ~new_rx
                end
                rs[a, 1] = new_rx
                if do_swap
                    rs[a, 2] = inv_ZR ? ~rx : rx
                elseif inv_ZR
                    rs[a, 2] = ~rz
                end
            else
                # X is multiplied into Z, so is the phase
                new_rz = rz ⊻ rx ⊻ pr
                if ~inv_ZR
                    new_rz = ~new_rz
                end
                rs[a, 2] = new_rz
                if do_swap
                    rs[a, 1] = inv_XR ? ~rz : rz
                elseif inv_XR
                    rs[a, 1] = ~rx
                end
            end
        end
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
        rs[a, 1] = u8_to_bool(rs[a, 1]) ⊻ u8_to_bool(rs[b, 1]) ⊻ (prod_phase_1 != 0)
        rs[b, 2] = u8_to_bool(rs[b, 2]) ⊻ u8_to_bool(rs[a, 2]) ⊻ (prod_phase_2 != 0)
    end
    return state
end

Base.@propagate_inbounds @inline function measure_z!(state::InvStabilizerState, a;
                                                     force=nothing)
    n = state.n
    @boundscheck check_qubit_bound(n, a)
    return _measure_z!(state, n, a, force)
end

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
    return @inbounds(u8_to_bool(rs[a, 2])), true

    @label rand_measure
    # Randomly pick a result
    res = force !== nothing ? force : rand(Bool)

    @inbounds begin
        xzs[lane + pchunk0, 3, a] = zero(VT2)
        cnot_tgt0 = zxa & ~pmask0
        za0 = xzs[lane + pchunk0, 4, a]
        xzs[lane + pchunk0, 4, a] = za0 | pmask0
        pza0 = _getbit(za0, pmask0)
        cnot_za0 = za0 & cnot_tgt0
    end

    i0_start = (pchunk0 >> 1) + 2
    if i0_start > i0_end
        @goto single_chunk
    end
    ws = state.ws
    assume(size(ws, 1) == 4)
    chunk_count = 0
    zcum_lo = cnot_za0
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
        ws[3, chunk_count] = i % ChT

        z = xzs[lane + i, 4, a]
        cnot_z = z & cnot_tgt
        zcum_hi ⊻= zcum_lo & cnot_z
        zcum_lo ⊻= cnot_z
    end
    if chunk_count <= 0
        @goto single_chunk
    end
    @inbounds begin
        zcum_lo_cnt_u8 = reduce(+, vcount_ones_u8(zcum_lo))
        was_y = pza0 ⊻ (zcum_lo_cnt_u8 & 1 != 0)

        zcum_hi_cnt_u8 = reduce(+, vcount_ones_u8(zcum_hi))
        zcum_cnt_u8 = (zcum_hi_cnt_u8 << 1) + zcum_lo_cnt_u8

        cnot_phase = zcum_lo_cnt_u8 ⊻ (ifelse(pza0, zcum_cnt_u8,
                                              zcum_cnt_u8 + 0x1) >> 1)
        flip_res = res ⊻ u8_to_bool(rs[a, 2]) ⊻ (cnot_phase & 1 != 0)
        rs[a, 2] = res
    end

    @inbounds for j in 1:n
        for k in 1:2
            if k == 2 && j == a
                continue
            end
            x0 = xzs[lane + pchunk0, 2k - 1, j]
            px = _getbit(x0, pmask0)
            z0 = xzs[lane + pchunk0, 2k, j]
            pz = _getbit(z0, pmask0)
            cnot_z0 = z0 & cnot_tgt0
            zcum_lo = cnot_z0
            if px
                zcum_hi = zero(VT2)
                xzcum = cnot_z0 & x0
                for cid in 1:chunk_count
                    cnot_tgt = ws[lane + 1, cid]
                    i = ws[3, cid] % Int
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
                r = u8_to_bool(rs[j, k]) ⊻ (cnot_phase & 1 != 0)
            else
                for cid in 1:chunk_count
                    cnot_tgt = ws[lane + 1, cid]
                    i = ws[3, cid] % Int
                    z = xzs[lane + i, 2k, j]
                    zcum_lo ⊻= z & cnot_tgt
                end
                zcum_lo_cnt_u8 = reduce(+, vcount_ones_u8(zcum_lo))
                r = u8_to_bool(rs[j, k])
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

    @inbounds if reduce(|, cnot_tgt0) != 0
        xzcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_za0 & zxa))
        zcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_za0))
        was_y = pza0 ⊻ (zcum_cnt_u8 & 1 != 0)
        cnot_phase = xzcum_cnt_u8 ⊻ (ifelse(pza0, zcum_cnt_u8,
                                            zcum_cnt_u8 + 0x1) >> 1)
        flip_res = res ⊻ u8_to_bool(rs[a, 2]) ⊻ (cnot_phase & 1 != 0)
        rs[a, 2] = res
        for j in 1:n
            for k in 1:2
                if k == 2 && j == a
                    continue
                end
                x = xzs[lane + pchunk0, 2k - 1, j]
                z = xzs[lane + pchunk0, 2k, j]
                px = _getbit(x, pmask0)
                pz = _getbit(z, pmask0)

                cnot_z = z & cnot_tgt0
                zcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_z))
                new_pz = pz ⊻ (zcum_cnt_u8 & 1 != 0)
                final_pz = ifelse(was_y, px, px ⊻ new_pz)
                final_px = ifelse(was_y, px ⊻ new_pz, new_pz)
                r = u8_to_bool(rs[j, k]) ⊻ (final_pz & flip_res)
                if px
                    xzcum_cnt_u8 = reduce(+, vcount_ones_u8(cnot_z & x))
                    cnot_phase = xzcum_cnt_u8 ⊻ (ifelse(pz, zcum_cnt_u8,
                                                        zcum_cnt_u8 + 0x1) >> 1)
                    x = x ⊻ cnot_tgt0
                    r ⊻= cnot_phase & 1 != 0
                end
                xzs[lane + pchunk0, 2k - 1, j] = _setbit(x, final_px, pmask0)
                xzs[lane + pchunk0, 2k, j] = _setbit(z, final_pz, pmask0)
                rs[j, k] = r
            end
        end
    else
        was_y = pza0
        flip_res = res ⊻ u8_to_bool(rs[a, 2])
        rs[a, 2] = res

        for j in 1:n
            for k in 1:2
                if k == 2 && j == a
                    continue
                end
                x = xzs[lane + pchunk0, 2k - 1, j]
                z = xzs[lane + pchunk0, 2k, j]

                px = x & pmask0
                pz = z & pmask0
                final_pz = ifelse(was_y, px, px ⊻ pz)
                flip_px = ifelse(was_y, pz, pz ⊻ px)
                xzs[lane + pchunk0, 2k - 1, j] = x ⊻ flip_px
                xzs[lane + pchunk0, 2k, j] = z ⊻ pz ⊻ final_pz
                rs[j, k] = (u8_to_bool(rs[j, k]) ⊻
                    ((reduce(|, final_pz) != 0) & flip_res))
            end
        end
    end
    return res, false
end

const _StabilizerState = Union{StabilizerState,InvStabilizerState}

@inline init_state_z!(state::_StabilizerState) = init_state_z!(state, false)
@inline init_state_x!(state::_StabilizerState) = init_state_x!(state, false)
@inline init_state_y!(state::_StabilizerState) = init_state_y!(state, false)

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

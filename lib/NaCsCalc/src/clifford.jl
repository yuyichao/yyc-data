#!/usr/bin/julia

module Clifford

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
Base.@propagate_inbounds function apply!(::CNOTGate, xas, zas, xbs, zbs, rs)
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

function Base.empty!(str::PauliString{T}) where T
    str.xs .= zero(T)
    str.zs .= zero(T)
    str.rs[] = zero(T)
    return str
end

function apply!(str::PauliString, gate::Clifford1Q, a)
    apply!(gate, @view(str.xs[a]), @view(str.zs[a]), str.rs)
    return str
end

function apply!(str::PauliString, gate::Clifford2Q, a, b)
    apply!(gate, @view(str.xs[a]), @view(str.zs[a]),
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
Base.@propagate_inbounds @inline function pauli_inject!(str::PauliString, x, z, i)
    str.xs[i] ⊻= x
    str.zs[i] ⊻= z
    return str
end
Base.@propagate_inbounds @inline inject_error!(str::PauliString, x, z, i) =
    pauli_inject!(str, x, z, i)

struct StabilizerState
    n::Int
    xs::Vector{BitVector}
    zs::Vector{BitVector}
    rs::BitVector
    function StabilizerState(n)
        xs = [(x = falses(2 * n + 1); x[i] = true; x) for i in 1:n]
        zs = [(z = falses(2 * n + 1); z[i + n] = true; z) for i in 1:n]
        rs = falses(2 * n + 1)
        return new(n, xs, zs, rs)
    end
end

function Base.show(io::IO, state::StabilizerState)
    for i in 1:state.n
        show(io, get_stabilizer(state, i))
        println(io)
    end
end

function get_stabilizer(state::StabilizerState, i)
    n = state.n
    return PauliString(n, [state.xs[j][i + n] for j in 1:n],
                       [state.zs[j][i + n] for j in 1:n], state.rs[i + n])
end

function init_state_0!(state::StabilizerState)
    for (i, x) in enumerate(state.xs)
        x.chunks .= 0
        x[i] = true
    end
    for (i, z) in enumerate(state.zs)
        z.chunks .= 0
        z[i + state.n] = true
    end
    state.rs.chunks .= 0
    return state
end

@inline function pauli_prod_phase(x1::Bool, z1::Bool, x2::Bool, z2::Bool)
    return ifelse(x1, ifelse(z1, z2 - x2, z2 * (2 * x2 - 1)),
                  ifelse(z1, x2 * (1 - 2 * z2), 0))
end
function clifford_rowsum!(state::StabilizerState, h, i)
    gs = 0
    for j in 1:state.n
        gs += pauli_prod_phase(state.xs[j][i], state.zs[j][i],
                               state.xs[j][h], state.zs[j][h])
    end
    gs += 2 * (state.rs[h] + state.rs[i])
    state.rs[h] = (gs & 3) != 0
    for j in 1:state.n
        state.xs[j][h] = state.xs[j][i] ⊻ state.xs[j][h]
        state.zs[j][h] = state.zs[j][i] ⊻ state.zs[j][h]
    end
    return state
end

# Randomly pick a result
function measure_z!(state::StabilizerState, a; force=nothing)
    n = state.n
    local p
    found_p = false
    for _p in (n + 1):(2 * n)
        if state.xs[a][_p]
            p = _p
            found_p = true
            break
        end
    end
    if found_p
        for i in 1:(2 * n)
            if i != p && state.xs[a][i]
                clifford_rowsum!(state, i, p)
            end
        end
        for i in 1:n
            state.xs[i][p - n] = state.xs[i][p]
            state.xs[i][p] = false
            state.zs[i][p - n] = state.zs[i][p]
            state.zs[i][p] = i == a
        end
        res = force !== nothing ? force : rand(Bool)
        state.rs[p - n] = state.rs[p]
        state.rs[p] = res
        return res, false
    else
        for i in 1:n
            state.xs[i][2 * n + 1] = false
            state.zs[i][2 * n + 1] = false
        end
        state.rs[2 * n + 1] = false
        for i in 1:n
            if state.xs[a][i]
                clifford_rowsum!(state, 2 * n + 1, i + n)
            end
        end
        return state.rs[2 * n + 1], true
    end
end

function measure_x!(state::StabilizerState, a; force=nothing)
    apply!(state, HGate(), a)
    res = measure_z!(state, a; force=force)
    apply!(state, HGate(), a)
    return res
end

function measure_y!(state::StabilizerState, a; force=nothing)
    apply!(state, SXGate(), a)
    res = measure_z!(state, a; force=force)
    apply!(state, ISXGate(), a)
    return res
end

function measure_zs!(state::StabilizerState, idxs; force=nothing)
    if isempty(idxs)
        return true, true
    end
    idx0 = idxs[1]
    nidxs = length(idxs)
    if nidxs == 1
        return measure_z!(state, idx0; force=force)
    end
    for i in 2:nidxs
        apply!(state, CNOTGate(), idxs[i], idx0)
    end
    res = measure_z!(state, idx0; force=force)
    for i in 2:nidxs
        apply!(state, CNOTGate(), idxs[i], idx0)
    end
    return res
end

function measure_xs!(state::StabilizerState, idxs; force=nothing)
    for idx in idxs
        apply!(state, HGate(), idx)
    end
    res = measure_zs!(state, idxs; force=force)
    for idx in idxs
        apply!(state, HGate(), idx)
    end
    return res
end

function measure_ys!(state::StabilizerState, idxs; force=nothing)
    for idx in idxs
        apply!(state, SXGate(), idx)
    end
    res = measure_zs!(state, idxs; force=force)
    for idx in idxs
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
        return true, true
    end
    res = measure_z!(state, idx0; force=force)
    for i in state.n:-1:1
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

function apply!(state::StabilizerState, gate::Clifford1Q, a)
    xas = state.xs[a]
    zas = state.zs[a]
    rs = state.rs
    nchunks = length(state.rs.chunks)
    @inbounds for i in 1:nchunks
        apply!(gate, @view(xas.chunks[i]), @view(zas.chunks[i]),
               @view(rs.chunks[i]))
    end
    return state
end

function apply!(state::StabilizerState, gate::Clifford2Q, a, b)
    xas = state.xs[a]
    zas = state.zs[a]
    xbs = state.xs[b]
    zbs = state.zs[b]
    rs = state.rs
    nchunks = length(state.rs.chunks)
    @inbounds for i in 1:nchunks
        apply!(gate, @view(xas.chunks[i]), @view(zas.chunks[i]),
               @view(xbs.chunks[i]), @view(zbs.chunks[i]),
               @view(rs.chunks[i]))
    end
    return state
end

@inline function inject_error!(state::StabilizerState, x::Bool, z::Bool, i)
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
    commute = zero(T)
    for (xa, za, xb, zb) in zip(xas, zas, xbs, zbs)
        xa = _cast_bits(T, xa)
        za = _cast_bits(T, za)
        xb = _cast_bits(T, xb)
        zb = _cast_bits(T, zb)
        commute ⊻= (xb & za) ⊻ (xa & zb)
    end
    return ~commute
end
pauli_commute(stra::PauliString, xbs, zbs) =
    pauli_commute(stra.xs, stra.zs, xbs, zbs)
pauli_commute(xas, zas, strb::PauliString) =
    pauli_commute(xas, zas, strb.xs, strb.zs)
pauli_commute(stra::PauliString, strb::PauliString) =
    pauli_commute(stra.xs, stra.zs, strb.xs, strb.zs)

end

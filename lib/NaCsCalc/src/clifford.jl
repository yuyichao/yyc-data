#!/usr/bin/julia

module Clifford

abstract type Clifford1Q end

abstract type Clifford2Q end

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

struct IGate <: Clifford1Q
end
@inline function apply!(::IGate, xas, zas, rs)
    return
end
@inline Base.inv(::IGate) = IGate()
Base.:*(g1::Clifford1Q, ::IGate) = g1
Base.:*(::IGate, g2::Clifford1Q) = g2
Base.:*(::IGate, ::IGate) = IGate()

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

struct PauliString
    n::Int
    xs::Vector{Bool}
    zs::Vector{Bool}
    rs::Base.RefValue{Bool}
    function PauliString(n)
        xs = fill(false, n)
        zs = fill(false, n)
        rs = Ref(false)
        return new(n, xs, zs, rs)
    end
    function PauliString(n, xs, zs, rs)
        return new(n, xs, zs, rs)
    end
end

function Base.empty!(state::PauliString)
    state.xs .= false
    state.zs .= false
    state.rs[] = false
    return state
end

function apply!(state::PauliString, gate::Clifford1Q, a)
    apply!(gate, @view(state.xs[a]), @view(state.zs[a]), state.rs)
    return state
end

function apply!(state::PauliString, gate::Clifford2Q, a, b)
    apply!(gate, @view(state.xs[a]), @view(state.zs[a]),
           @view(state.xs[b]), @view(state.zs[b]), state.rs)
    return state
end

function Base.show(io::IO, str::PauliString)
    write(io, str.rs[] ? '-' : '+')
    for (x, z) in zip(str.xs, str.zs)
        write(io, x ? (z ? "Y" : "X") : (z ? "Z" : "I"))
    end
end

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

function get_stabilizer(state::StabilizerState, i)
    n = state.n
    return PauliString(n, [state.xs[j][i] for j in 1:n],
                       [state.zs[j][i] for j in 1:n], state.rs[i])
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

function pauli_commute(xas, zas, xbs, zbs)
    commute = zero(eltype(xas))
    for (xa, za, xb, zb) in zip(xas, zas, xbs, zbs)
        commute ⊻= (xb & za) ⊻ (xa & zb)
    end
    return ~commute
end

end

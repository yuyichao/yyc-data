#!/usr/bin/julia

module Clifford

abstract type Clifford1Q end

abstract type Clifford2Q end

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

struct XGate <: Clifford1Q
end
# I -> I, X -> X, Y -> -Y, Z -> -Z
Base.@propagate_inbounds @inline function apply!(::XGate, xas, zas, rs)
    za = zas[]
    r = rs[]

    rs[] = r ⊻ za
    return
end

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

struct ZGate <: Clifford1Q
end
# I -> I, X -> -X, Y -> -Y, Z -> Z
Base.@propagate_inbounds @inline function apply!(::ZGate, xas, zas, rs)
    xa = xas[]
    r = rs[]

    rs[] = r ⊻ xa
    return
end

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

end

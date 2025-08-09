#!/usr/bin/julia

mutable struct Leaf
    const resistance::Float64
    flow::Float64
    Leaf(r) = new(r)
end

mutable struct Branch
    const b1::Union{Branch,Leaf}
    const b2::Union{Branch,Leaf}
    const resistance0::Float64
    resistance::Float64
    flow::Float64
    Branch(b1, b2, r0) = new(b1, b2, r0)
end

function update_resistance!(b::Branch)
    r1 = update_resistance!(b.b1)
    r2 = update_resistance!(b.b2)
    b.resistance = r1 * r2 / (r1 + r2) + b.resistance0
    return b.resistance
end

function update_resistance!(l::Leaf)
    return l.resistance
end

function update_flow!(b::Branch, flow)
    r1 = b.b1.resistance
    r2 = b.b2.resistance
    b.flow = flow
    update_flow!(b.b1, flow * r2 / (r1 + r2))
    update_flow!(b.b2, flow * r1 / (r1 + r2))
end

function update_flow!(l::Leaf, flow)
    l.flow = flow
    return
end

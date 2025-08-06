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

const out_resistance = 100.0

const l1 = Leaf(5 + out_resistance)
const l2 = Leaf(6 + out_resistance)
const l3 = Leaf(2.5 + out_resistance)
const l4 = Leaf(4 + out_resistance)
const l5 = Leaf(10 + out_resistance)
const l6 = Leaf(1.5 + out_resistance)
const l7 = Leaf(2 + out_resistance)
const l8 = Leaf(5 + out_resistance)
const l9 = Leaf(2 + out_resistance)
const l10 = Leaf(9 + out_resistance)
const l11 = Leaf(4.5 + out_resistance)
const l12 = Leaf(2 + out_resistance)
const l13 = Leaf(4 + out_resistance)
const l14 = Leaf(1 + out_resistance)
const l15 = Leaf(8 + out_resistance)
const l16 = Leaf(3 + out_resistance)
const l17 = Leaf(0.5 + out_resistance)

const b1 = Branch(l3, l4, 18)
const b2 = Branch(l2, b1, 4)
const b3 = Branch(l1, b2, 5)
const b4 = Branch(b3, l5, 5)
const b5 = Branch(l9, l10, 8)
const b6 = Branch(l8, b5, 2)
const b7 = Branch(b6, l11, 13)
const b8 = Branch(l7, b7, 7)
const b9 = Branch(l6, b8, 6)
const b10 = Branch(b4, b9, 10)
const b11 = Branch(l12, l13, 7)
const b12 = Branch(l14, l15, 3)
const b13 = Branch(b12, l16, 4)
const b14 = Branch(b13, l17, 5)
const b15 = Branch(b11, b14, 37)
const b16 = Branch(b10, b15, 0)

const root = b16

update_resistance!(root)
update_flow!(root, 1)

println("L1: $(l1.flow * 100) %")
println("L2: $(l2.flow * 100) %")
println("L3: $(l3.flow * 100) %")
println("L4: $(l4.flow * 100) %")
println("L5: $(l5.flow * 100) %")
println("L6: $(l6.flow * 100) %")
println("L7: $(l7.flow * 100) %")
println("L8: $(l8.flow * 100) %")
println("L9: $(l9.flow * 100) %")
println("L10: $(l10.flow * 100) %")
println("L11: $(l11.flow * 100) %")
println("L12: $(l12.flow * 100) %")
println("L13: $(l13.flow * 100) %")
println("L14: $(l14.flow * 100) %")
println("L15: $(l15.flow * 100) %")
println("L16: $(l16.flow * 100) %")
println("L17: $(l17.flow * 100) %")

println("B1: $(b1.flow * 100) %")
println("B2: $(b2.flow * 100) %")
println("B3: $(b3.flow * 100) %")
println("B4: $(b4.flow * 100) %")
println("B5: $(b5.flow * 100) %")
println("B6: $(b6.flow * 100) %")
println("B7: $(b7.flow * 100) %")
println("B8: $(b8.flow * 100) %")
println("B9: $(b9.flow * 100) %")
println("B10: $(b10.flow * 100) %")
println("B11: $(b11.flow * 100) %")
println("B12: $(b12.flow * 100) %")
println("B13: $(b13.flow * 100) %")
println("B14: $(b14.flow * 100) %")
println("B15: $(b15.flow * 100) %")
println("B16: $(b16.flow * 100) %")

println("$(l1.flow * 100)")
println("$(l2.flow * 100)")
println("$(l3.flow * 100)")
println("$(l4.flow * 100)")
println("$(l5.flow * 100)")
println("$(l6.flow * 100)")
println("$(l7.flow * 100)")
println("$(l8.flow * 100)")
println("$(l9.flow * 100)")
println("$(l10.flow * 100)")
println("$(l11.flow * 100)")
println("$(l12.flow * 100)")
println("$(l13.flow * 100)")
println("$(l14.flow * 100)")
println("$(l15.flow * 100)")
println("$(l16.flow * 100)")
println("$(l17.flow * 100)")

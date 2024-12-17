#!/usr/bin/julia

const va0 = 30899381

@inline function eval_round(va)
    vb = va & 7
    vb = vb ⊻ 1
    vc = va >> vb
    vb = vb ⊻ vc
    va = va >> 3
    vb = vb ⊻ 6
    return vb & 7, va
end

function eval_all(va)
    outputs = Int[]
    while va != 0
        o, va = eval_round(va)
        push!(outputs, o)
    end
    return outputs
end

@show join(eval_all(va0), ",")

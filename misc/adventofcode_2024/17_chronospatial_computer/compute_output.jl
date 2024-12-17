#!/usr/bin/julia

const va0 = 30899381
const vb0 = 0
const vc0 = 0

const code = [2,4,1,1,7,5,4,0,0,3,1,6,5,5,3,0]

mutable struct Processor
    va::Int
    vb::Int
    vc::Int
    ip::Int
end

function get_combo_operand(p::Processor, i)
    if i <= 3
        return i
    elseif i == 4
        return p.va
    elseif i == 5
        return p.vb
    elseif i == 6
        return p.vc
    else
        error("...")
    end
end

function eval_inst(p::Processor, code, outputs)
    ip = p.ip
    if ip >= length(code)
        return false
    end
    opcode = code[ip + 1]
    operand = code[ip + 2]
    p.ip = ip + 2
    if opcode == 0
        p.va = p.va >> get_combo_operand(p, operand)
    elseif opcode == 1
        p.vb = p.vb ⊻ operand
    elseif opcode == 2
        p.vb = get_combo_operand(p, operand) & 7
    elseif opcode == 3
        if p.va != 0
            p.ip = operand
        end
    elseif opcode == 4
        p.vb = p.vb ⊻ p.vc
    elseif opcode == 5
        push!(outputs, get_combo_operand(p, operand) & 7)
    elseif opcode == 6
        p.vb = p.va >> get_combo_operand(p, operand)
    elseif opcode == 7
        p.vc = p.va >> get_combo_operand(p, operand)
    end
    return true
end

const outputs = Int[]
const p = Processor(va0, vb0, vc0, 0)
while eval_inst(p, code, outputs)
end
@show join(outputs, ",")

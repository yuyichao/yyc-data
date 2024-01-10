#!/usr/bin/julia

mutable struct PState
    accum::Int
    pc::Int
    PState() = new(0, 1)
end

@enum OpCode Nop Acc Jmp

struct Inst
    op::OpCode
    arg::Int
end

function load_insts(file)
    insts = Inst[]
    for line in eachline(file)
        m = match(r"(nop|acc|jmp) ([+-]\d+)", line)
        arg = parse(Int, m[2])
        if m[1] == "nop"
            op = Nop
        elseif m[1] == "acc"
            op = Acc
        elseif m[1] == "jmp"
            op = Jmp
        else
            error()
        end
        push!(insts, Inst(op, arg))
    end
    return insts
end

function find_loop(file)
    insts = load_insts(file)
    state = PState()
    seen_pc = Set{Int}()
    while true
        if state.pc in seen_pc
            return state.accum
        end
        push!(seen_pc, state.pc)
        inst = insts[state.pc]
        if inst.op == Nop
            state.pc += 1
        elseif inst.op == Acc
            state.pc += 1
            state.accum += inst.arg
        elseif inst.op == Jmp
            state.pc += inst.arg
        else
            error()
        end
    end
end

@show find_loop(ARGS[1])

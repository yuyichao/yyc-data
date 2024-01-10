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

function check_halt(insts)
    state = PState()
    ninsts = length(insts)

    for i in 1:ninsts
        if state.pc == ninsts + 1
            return state.accum
        end
        if !checkbounds(Bool, insts, state.pc)
            return
        end
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

function find_exit(file)
    insts = load_insts(file)

    ninsts = length(insts)
    for i in 1:ninsts
        inst = insts[i]
        if inst.op == Nop
            insts[i] = Inst(Jmp, inst.arg)
        elseif inst.op == Jmp
            insts[i] = Inst(Nop, inst.arg)
        else
            continue
        end
        v = check_halt(insts)
        if v !== nothing
            return v
        end
        insts[i] = inst
    end
end

@show find_exit(ARGS[1])

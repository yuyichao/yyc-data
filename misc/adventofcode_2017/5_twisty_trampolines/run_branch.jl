#!/usr/bin/julia

function run_inst(insts, idx)
    inst = insts[idx]
    insts[idx] = inst + 1
    idx += inst
    if idx < 1 || idx > length(insts)
        return
    end
    return idx
end

function find_escape(file)
    insts = [parse(Int, line) for line in eachline(file)]

    idx = 1
    c = 0
    while idx !== nothing
        c += 1
        idx = run_inst(insts, idx)
    end
    return c
end

@show find_escape(ARGS[1])

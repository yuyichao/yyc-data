#!/usr/bin/julia

using OrderedCollections

function hash_func(s)
    h = 0x0
    for c in s
        h += c % UInt8
        h *= UInt8(17)
    end
    return Int(h)
end

function read_instructions(file)
    s = read(file, String)
    s = strip(s)
    return split(s, ',')
end

function run_instructions(ins)
    state = [OrderedDict{String,Int}() for i in 1:256]
    for s in ins
        m = match(r"^([a-z]+)-$", s)
        if m !== nothing
            label = m[1]
            delete!(state[hash_func(label) + 1], label)
            continue
        end
        m = match(r"^([a-z]+)=([0-9]+)$", s)
        if m !== nothing
            label = m[1]
            state[hash_func(label) + 1][label] = parse(Int, m[2])
            continue
        end
        error()
    end
    return state
end

function total_focus(file)
    state = run_instructions(read_instructions(file))
    total_f = 0
    for i in 1:256
        lens_order = 1
        for f in values(state[i])
            total_f += i * lens_order * f
            lens_order += 1
        end
    end
    return total_f
end

@show total_focus(ARGS[1])

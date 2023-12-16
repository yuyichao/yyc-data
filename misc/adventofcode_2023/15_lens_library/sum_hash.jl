#!/usr/bin/julia

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

function sum_hash(file)
    ins = read_instructions(file)
    return sum(hash_func(s) for s in ins)
end

@show sum_hash(ARGS[1])

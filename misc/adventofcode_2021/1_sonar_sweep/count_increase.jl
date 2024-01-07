#!/usr/bin/julia

function count_increase(file)
    prev = typemax(Int)
    c = 0
    for line in eachline(file)
        v = parse(Int, line)
        c += v > prev
        prev = v
    end
    return c
end

@show count_increase(ARGS[1])

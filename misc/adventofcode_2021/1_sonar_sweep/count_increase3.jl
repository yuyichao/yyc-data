#!/usr/bin/julia

function count_increase(file)
    prev1, prev2, prev3 = typemax(Int), typemax(Int), typemax(Int)
    c = 0
    for line in eachline(file)
        v = parse(Int, line)
        c += v > prev1
        prev1, prev2, prev3 = prev2, prev3, v
    end
    return c
end

@show count_increase(ARGS[1])

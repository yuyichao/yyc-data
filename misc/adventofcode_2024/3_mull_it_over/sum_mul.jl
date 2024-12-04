#!/usr/bin/julia

function sum_mul(file)
    s = 0
    for line in eachline(file)
        for m in eachmatch(r"mul\(([0-9]+),([0-9]+)\)", line)
            s += parse(Int, m[1]) * parse(Int, m[2])
        end
    end
    return s
end

@show sum_mul(ARGS[1])

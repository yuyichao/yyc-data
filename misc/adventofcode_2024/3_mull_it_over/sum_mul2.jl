#!/usr/bin/julia

function sum_mul(file)
    s = 0
    enabled = true
    for line in eachline(file)
        for m in eachmatch(r"mul\(([0-9]+),([0-9]+)\)|do\(\)|don't\(\)", line)
            if m.match == "do()"
                enabled = true
                continue
            elseif m.match == "don't()"
                enabled = false
                continue
            end
            if enabled
                s += parse(Int, m[1]) * parse(Int, m[2])
            end
        end
    end
    return s
end

@show sum_mul(ARGS[1])

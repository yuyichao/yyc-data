#!/usr/bin/julia

function sum_file(file)
    s = 0
    for line in eachline(file)
        matches = collect(eachmatch(r"[0-9]", line))
        d1 = parse(Int, matches[1].match)
        d2 = parse(Int, matches[end].match)
        s += d1 * 10 + d2
    end
    return s
end

@show sum_file(ARGS[1])

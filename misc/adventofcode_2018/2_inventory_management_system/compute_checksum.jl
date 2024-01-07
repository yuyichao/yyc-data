#!/usr/bin/julia

using DataStructures

function classify_id(id)
    is2 = false
    is3 = false
    for (c, n) in counter(id)
        if n == 2
            is2 = true
        elseif n == 3
            is3 = true
        end
    end
    return is2, is3
end

function compute_checksum(file)
    c2 = 0
    c3 = 0
    for line in eachline(file)
        is2, is3 = classify_id(line)
        c2 += is2
        c3 += is3
    end
    return c2 * c3
end

@show compute_checksum(ARGS[1])

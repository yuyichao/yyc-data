#!/usr/bin/julia

function sum_part(file)
    numbers = Vector{NTuple{3,Int}}[] # value, start, end
    symbols = Vector{Bool}[]
    for line in eachline(file)
        push!(numbers, [(parse(Int, m.match), m.offset - 1,
                         m.offset + length(m.match))
                        for m in eachmatch(r"[0-9]+", line)])
        push!(symbols, [(s < '0' || s > '9') && s != '.' for s in line])
    end

    nlines = length(numbers)
    ncols = length(symbols[1])
    s = 0
    for lo in 1:nlines
        for (num, i1, i2) in numbers[lo]
            i1 = max(i1, 1)
            i2 = min(i2, ncols)
            for l in lo - 1:lo + 1
                if l < 1 || l > nlines
                    continue
                end
                if any(symbols[l][i1:i2])
                    s += num
                    break
                end
            end
        end
    end
    return s
end

@show sum_part(ARGS[1])

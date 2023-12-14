#!/usr/bin/julia

function sum_gear(file)
    numbers = Vector{NTuple{3,Int}}[] # value, start, end
    gears = Vector{Int}[]
    for line in eachline(file)
        push!(numbers, [(parse(Int, m.match), m.offset - 1,
                         m.offset + length(m.match))
                        for m in eachmatch(r"[0-9]+", line)])
        push!(gears, first.(findall("*", line)))
    end

    nlines = length(numbers)
    s = 0
    for lo in 1:nlines
        for gi in gears[lo]
            num_count = 0
            ratio = 1
            for l in lo - 1:lo + 1
                if l < 1 || l > nlines
                    continue
                end
                for (num, i1, i2) in numbers[l]
                    if i1 <= gi <= i2
                        ratio *= num
                        num_count += 1
                    end
                    if num_count > 2
                        break
                    end
                end
                if num_count > 2
                    break
                end
            end
            if num_count == 2
                s += ratio
            end
        end
    end
    return s
end

@show sum_gear(ARGS[1])

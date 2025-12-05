#!/usr/bin/julia

function count_fresh(file)
    lines = eachline(file)
    ranges = Set{UnitRange{Int}}()
    for line in lines
        if isempty(line)
            break
        end
        lb, ub = parse.(Int, split(line, "-"))
        push!(ranges, lb:ub)
    end
    c = 0
    for line in lines
        v = parse(Int, line)
        for rng in ranges
            if v in rng
                c += 1
                break
            end
        end
    end
    return c
end

@show count_fresh(ARGS[1])

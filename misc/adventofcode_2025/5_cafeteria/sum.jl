#!/usr/bin/julia

function count_fresh(file)
    lines = eachline(file)
    ranges = NTuple{2,Int}[]
    for line in lines
        if isempty(line)
            break
        end
        lb, ub = parse.(Int, split(line, "-"))
        push!(ranges, (lb, ub))
    end
    sort!(ranges)
    cur_lb, cur_ub = ranges[1]
    c = 0
    for i in 2:length(ranges)
        lb, ub = ranges[i]
        if lb > cur_ub + 1
            c += cur_ub - cur_lb + 1
            cur_lb, cur_ub = lb, ub
        else
            cur_ub = max(ub, cur_ub)
        end
    end
    c += cur_ub - cur_lb + 1
    return c
end

@show count_fresh(ARGS[1])

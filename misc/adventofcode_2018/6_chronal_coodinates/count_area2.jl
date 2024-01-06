#!/usr/bin/julia

function count_area(file)
    points = NTuple{2,Int}[]
    min_x = min_y = typemax(Int)
    max_x = max_y = typemin(Int)
    for line in eachline(file)
        x, y = parse.(Int, split(line, ','))
        min_x = min(x, min_x)
        min_y = min(y, min_y)
        max_x = max(x, max_x)
        max_y = max(y, max_y)
        push!(points, (x, y))
    end
    c = 0
    for x in min_x:max_x
        for y in min_y:max_y
            d = sum(p->abs(x - p[1]) + abs(y - p[2]), points)
            if d < 10000
                c += 1
            end
        end
    end
    return c
end

@show count_area(ARGS[1])

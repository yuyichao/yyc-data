#!/usr/bin/julia

function find_closest_point(points, x, y)
    min_d = typemax(Int)
    min_p = (typemax(Int), typemax(Int))
    for p in points
        d = abs(x - p[1]) + abs(y - p[2])
        if min_d > d
            min_d = d
            min_p = p
        elseif min_d == d
            min_p = (typemax(Int), typemax(Int))
        end
    end
    return min_p
end

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
    counts = Dict{NTuple{2,Int},Int}()
    for x in min_x:max_x
        for y in min_y:max_y
            p = find_closest_point(points, x, y)
            counts[p] = get(counts, p, 0) + 1
        end
    end

    delete!(counts, (typemax(Int), typemax(Int)))
    for x in (min_x - 1):(max_x + 1)
        delete!(counts, find_closest_point(points, x, min_y - 1))
        delete!(counts, find_closest_point(points, x, max_y + 1))
    end
    for y in (min_y - 1):(max_y + 1)
        delete!(counts, find_closest_point(points, min_x - 1, y))
        delete!(counts, find_closest_point(points, max_x + 1, y))
    end

    return maximum(x->x[2], counts)
end

@show count_area(ARGS[1])

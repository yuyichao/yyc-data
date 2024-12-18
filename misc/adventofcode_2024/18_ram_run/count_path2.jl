#!/usr/bin/julia

function check_coord(M, x, y, newd, todo)
    if !checkbounds(Bool, M, x, y)
        return
    end
    prev = M[x, y]
    if prev == -1 || prev > newd
        M[x, y] = newd
        push!(todo, (x, y))
    end
    return
end

function count_path(coords)
    M = fill(-1, 71, 71)
    for (x, y) in coords
        M[x + 1, y + 1] = -2
    end
    M[1, 1] = 0
    todo = Set(((1, 1),))
    while !isempty(todo)
        x, y = pop!(todo)
        d = M[x, y]
        newd = d + 1
        check_coord(M, x + 1, y, newd, todo)
        check_coord(M, x - 1, y, newd, todo)
        check_coord(M, x, y + 1, newd, todo)
        check_coord(M, x, y - 1, newd, todo)
    end
    return M[71, 71]
end

function min_drop(file)
    coords = NTuple{2,Int}[]
    for line in eachline(file)
        x, y = parse.(Int, split(line, ","))
        push!(coords, (x, y))
    end
    @assert count_path(coords) < 0
    lo = 1024
    hi = length(coords)
    while hi > lo + 1
        mid = (lo + hi) ÷ 2
        if count_path(coords[1:mid]) < 0
            hi = mid
        else
            lo = mid
        end
    end
    return coords[hi]
end

@show min_drop(ARGS[1])

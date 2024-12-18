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

function count_path(file)
    lines = readlines(file)[1:1024]
    M = fill(-1, 71, 71)
    for line in lines
        x, y = parse.(Int, split(line, ","))
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

@show count_path(ARGS[1])

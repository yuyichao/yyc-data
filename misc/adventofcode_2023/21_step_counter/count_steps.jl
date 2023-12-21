#!/usr/bin/julia

function propagate!(garden_map, x0, y0, max_step, frontier)
    step = garden_map[y0, x0]
    next_step = step + 1
    if next_step > max_step
        return
    end
    for (x, y) in ((x0 + 1, y0), (x0 - 1, y0), (x0, y0 + 1), (x0, y0 - 1))
        if !checkbounds(Bool, garden_map, y, x) || garden_map[y, x] <= next_step
            continue
        end
        garden_map[y, x] = next_step
        push!(frontier, (x, y))
    end
end

function propagate_all!(garden_map, max_step, frontier)
    while !isempty(frontier)
        (x, y) = popfirst!(frontier)
        propagate!(garden_map, x, y, max_step, frontier)
    end
end

function count_steps(file)
    frontier = NTuple{2,Int}[]

    lines = readlines(file)
    M = Matrix{Int}(undef, length(lines), length(lines[1]))
    for i in 1:length(lines)
        line = lines[i]
        for j in 1:length(line)
            c = line[j]
            if c == '#'
                M[i, j] = -1
            elseif c == 'S'
                M[i, j] = 0
                push!(frontier, (j, i))
            else
                @assert c == '.'
                M[i, j] = typemax(Int)
            end
        end
    end
    propagate_all!(M, 64, frontier)
    return count(x->x <= 64 && x % 2 == 0, M)
end

@show count_steps(ARGS[1])

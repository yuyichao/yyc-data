#!/usr/bin/julia

function load_map(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int8}(undef, nrow, ncol)
    start = (0, 0)
    for i in 1:nrow
        line = lines[i]
        for j in 1:ncol
            c = line[j]
            if c == '.'
                M[i, j] = 0
            elseif c == '^'
                start = (i, j)
                M[i, j] = 0
            elseif c == '#'
                M[i, j] = 1
            else
                error("Unknown character")
            end
        end
    end
    return M, start
end

function count_area(file)
    M, pos = load_map(file)
    dir = (-1, 0)
    c = 0
    while true
        if M[pos...] == 0
            M[pos...] = -1
            c += 1
        end
        next_pos = pos .+ dir
        if !checkbounds(Bool, M, next_pos...)
            break
        end
        if M[next_pos...] != 1
            pos = next_pos
        else
            dir = (dir[2], -dir[1])
        end
    end
    return c
end

@show count_area(ARGS[1])

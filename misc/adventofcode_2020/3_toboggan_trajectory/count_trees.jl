#!/usr/bin/julia

function read_map(file)
    lines = readlines(file)
    nrows = length(lines)
    ncols = length(lines[1])
    M = zeros(Bool, nrows, ncols)
    for y in 1:nrows
        line = lines[y]
        for x in 1:ncols
            M[y, x] = line[x] == '#'
        end
    end
    return M
end

function count_trees_slope(M, slope)
    x = 1
    c = 0
    for y in 1:size(M, 1)
        if M[y, x]
            c += 1
        end
        x = (x + slope - 1) % size(M, 2) + 1
    end
    return c
end

@show count_trees_slope(read_map(ARGS[1]), 3)

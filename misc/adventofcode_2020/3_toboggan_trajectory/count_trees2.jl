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

function count_trees_step(M, xstep, ystep)
    x = 1
    y = 1
    nrows, ncols = size(M)

    c = 0
    while y <= nrows
        if M[y, x]
            c += 1
        end
        x = (x + xstep - 1) % size(M, 2) + 1
        y += ystep
    end

    return c
end

function count_trees(file)
    M = read_map(file)

    return count_trees_step(M, 1, 1) * count_trees_step(M, 3, 1) * count_trees_step(M, 5, 1) * count_trees_step(M, 7, 1) * count_trees_step(M, 1, 2)
end

@show count_trees(ARGS[1])

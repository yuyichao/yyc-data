#!/usr/bin/julia

function find_basin_from(M, i, j)
    points = Set(((i, j),))
    todo = copy(points)
    while !isempty(todo)
        p0 = pop!(todo)
        for p in (p0 .+ (1, 0), p0 .+ (-1, 0), p0 .+ (0, 1), p0 .+ (0, -1))
            if checkbounds(Bool, M, p...) && M[p...] != 9 && !(p in points)
                push!(todo, p)
                push!(points, p)
            end
        end
    end
    return length(points)
end

function find_basins(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int8}(undef, nrow, ncol)
    for i in 1:nrow
        M[i, :] .= (c - '0' for c in lines[i])
    end

    sizes = Int[]

    for i in 1:nrow
        for j in 1:ncol
            h = M[i, j]
            if i > 1 && h >= M[i - 1, j]
                continue
            elseif i < nrow && h >= M[i + 1, j]
                continue
            elseif j > 1 && h >= M[i, j - 1]
                continue
            elseif j < ncol && h >= M[i, j + 1]
                continue
            end
            push!(sizes, find_basin_from(M, i, j))
        end
    end

    sort!(sizes)

    # @show (sizes[end], sizes[end - 1], sizes[end - 2])
    return sizes[end] * sizes[end - 1] * sizes[end - 2]
end

@show find_basins(ARGS[1])

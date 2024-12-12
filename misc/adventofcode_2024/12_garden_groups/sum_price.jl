#!/usr/bin/julia

function area_price(M, assigned, row, col)
    c = M[row, col]
    if assigned[row, col]
        return 0
    end
    a = 1
    assigned[row, col] = true
    l = 0
    todo = Set(((row, col),))
    while !isempty(todo)
        p0 = pop!(todo)
        for p in (p0 .+ (1, 0), p0 .+ (-1, 0), p0 .+ (0, 1), p0 .+ (0, -1))
            if !checkbounds(Bool, M, p...) || M[p...] != c
                l += 1
            elseif !assigned[p...]
                a += 1
                assigned[p...] = true
                push!(todo, p)
            end
        end
    end
    return a * l
end

function sum_price(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int}(undef, nrow, ncol)
    for row in 1:nrow
        line = lines[row]
        M[row, :] .= (c for c in line)
    end
    assigned = zeros(Bool, nrow, ncol)
    s = 0
    for row in 1:nrow
        for col in 1:ncol
            s += area_price(M, assigned, row, col)
        end
    end
    return s
end

@show sum_price(ARGS[1])

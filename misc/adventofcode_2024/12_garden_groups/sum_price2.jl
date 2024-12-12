#!/usr/bin/julia

function count_sides(side)
    s = 0
    for coord in values(side)
        sort!(coord)
        last = -1
        for v in coord
            if v != last + 1
                s += 1
            end
            last = v
        end
    end
    return s
end

function area_price(M, assigned, row, col)
    c = M[row, col]
    if assigned[row, col]
        return 0
    end
    a = 1
    assigned[row, col] = true
    todo = Set(((row, col),))
    sides = (Dict{Int,Vector{Int}}(),
             Dict{Int,Vector{Int}}(),
             Dict{Int,Vector{Int}}(),
             Dict{Int,Vector{Int}}())
    while !isempty(todo)
        p0 = pop!(todo)
        ps = (p0 .+ (1, 0), p0 .+ (-1, 0), p0 .+ (0, 1), p0 .+ (0, -1))
        for i in 1:4
            p = ps[i]
            if !checkbounds(Bool, M, p...) || M[p...] != c
                if i == 1 || i == 2
                    push!(get!(()->Int[], sides[i], p0[1]), p0[2])
                else
                    push!(get!(()->Int[], sides[i], p0[2]), p0[1])
                end
            elseif !assigned[p...]
                a += 1
                assigned[p...] = true
                push!(todo, p)
            end
        end
    end
    return a * sum(count_sides(side) for side in sides)
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

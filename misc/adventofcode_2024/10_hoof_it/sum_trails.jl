#!/usr/bin/julia

function get_target(M, i, j, targets)
    tgt = targets[i, j]
    if !isempty(tgt)
        return tgt
    end
    v = M[i, j]
    if v == 9
        push!(tgt, (i, j))
        return tgt
    end
    for (i2, j2) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))
        if !checkbounds(Bool, M, i2, j2)
            continue
        end
        if M[i2, j2] != v + 1
            continue
        end
        union!(tgt, get_target(M, i2, j2, targets))
    end
    return tgt
end

function sum_trails(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int}(undef, nrow, ncol)
    for i in 1:nrow
        line = lines[i]
        M[i, :] .= (c - '0' for c in line)
    end
    targets = [Set{NTuple{2,Int}}() for i in 1:nrow, j in 1:ncol]
    s = 0
    for i in 1:nrow
        for j in 1:ncol
            if M[i, j] == 0
                s += length(get_target(M, i, j, targets))
            end
        end
    end
    return s
end

@show sum_trails(ARGS[1])

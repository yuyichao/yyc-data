#!/usr/bin/julia

function get_score(M, i, j, scores)
    if scores[i, j] != -1
        return scores[i, j]
    end
    v = M[i, j]
    if v == 9
        scores[i, j] = 1
        return 1
    end
    s = 0
    for (i2, j2) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))
        if !checkbounds(Bool, M, i2, j2)
            continue
        end
        if M[i2, j2] != v + 1
            continue
        end
        s += get_score(M, i2, j2, scores)
    end
    scores[i, j] = s
    return s
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
    scores = fill(-1, (nrow, ncol))
    s = 0
    for i in 1:nrow
        for j in 1:ncol
            if M[i, j] == 0
                s += get_score(M, i, j, scores)
            end
        end
    end
    return s
end

@show sum_trails(ARGS[1])

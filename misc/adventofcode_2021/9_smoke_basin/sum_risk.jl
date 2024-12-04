#!/usr/bin/julia

function sum_risk(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int8}(undef, nrow, ncol)
    for i in 1:nrow
        M[i, :] .= (c - '0' for c in lines[i])
    end
    s = 0
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
            s += h + 1
        end
    end
    return s
end

@show sum_risk(ARGS[1])

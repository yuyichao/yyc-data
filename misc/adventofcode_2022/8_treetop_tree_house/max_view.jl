#!/usr/bin/julia

function get_view(M, i0, j0)
    res = 1
    h0 = M[i0, j0]
    for delta in ((1, 0), (-1, 0), (0, 1), (0, -1))
        d = 0
        while true
            d += 1
            i, j = (i0, j0) .+ delta .* d
            if !checkbounds(Bool, M, i, j)
                res *= (d - 1)
                break
            elseif M[i, j] >= h0
                res *= d
                break
            end
        end
    end
    return res
end

function max_view(file)
    lines = readlines(file)

    nrows = length(lines)
    ncols = length(lines[1])

    M = Matrix{Int}(undef, ncols, nrows)

    for i in 1:nrows
        line = lines[i]
        for j in 1:ncols
            M[i, j] = line[j] - '0'
        end
    end

    return maximum(get_view(M, i, j) for i in 1:nrows, j in 1:ncols)
end

@show max_view(ARGS[1])

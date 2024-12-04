#!/usr/bin/julia

function get_idx(M, args...)
    if !checkbounds(Bool, M, args...)
        return zero(eltype(M))
    end
    return M[args...]
end

function count_word(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int8}(undef, nrow, ncol)
    for i in 1:nrow
        M[i, :] .= (c for c in lines[i])
    end

    target = Int8.(('M', 'A', 'S'))
    function get_didx(i, j, di, dj)
        return (get_idx(M, i + di, j + dj),
                get_idx(M, i + di * 2, j + dj * 2),
                get_idx(M, i + di * 3, j + dj * 3))
    end

    s = 0
    for i in 1:nrow
        for j in 1:ncol
            h = M[i, j]
            if h != Int8('X')
                continue
            end
            for di in -1:1
                for dj in -1:1
                    if di == dj == 0
                        continue
                    end
                    if get_didx(i, j, di, dj) == target
                        s += 1
                    end
                end
            end
        end
    end
    return s
end

@show count_word(ARGS[1])

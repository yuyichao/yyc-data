#!/usr/bin/julia

function get_idx(M, args...)
    if !checkbounds(Bool, M, args...)
        return zero(eltype(M))
    end
    return M[args...]
end

function count_x(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int8}(undef, nrow, ncol)
    for i in 1:nrow
        M[i, :] .= (c for c in lines[i])
    end

    # target = Int8.(('M', 'A', 'S'))
    # function get_didx(i, j, di, dj)
    #     return (get_idx(M, i + di, j + dj),
    #             get_idx(M, i + di * 2, j + dj * 2),
    #             get_idx(M, i + di * 3, j + dj * 3))
    # end

    targets = (Int8('M'), Int8('S')), (Int8('S'), Int8('M'))

    s = 0
    for i in 2:(nrow - 1)
        for j in 2:(ncol - 1)
            h = M[i, j]
            if h != Int8('A')
                continue
            end
            d1 = M[i - 1, j - 1], M[i + 1, j + 1]
            d2 = M[i - 1, j + 1], M[i + 1, j - 1]

            if d1 in targets && d2 in targets
                s += 1
            end
        end
    end
    return s
end

@show count_x(ARGS[1])

#!/usr/bin/julia

function load_map(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int8}(undef, nrow, ncol)
    start = (0, 0)
    for i in 1:nrow
        line = lines[i]
        for j in 1:ncol
            c = line[j]
            if c == '.'
                M[i, j] = 0
            elseif c == '^'
                start = (i, j)
                M[i, j] = 0
            elseif c == '#'
                M[i, j] = -1
            else
                error("Unknown character")
            end
        end
    end
    return M, start
end

function check_loop(M, pos, dir)
    dir_dpos = ((-1, 0), (0, 1), (1, 0), (0, -1))
    while true
        mask = 1 << dir
        if (M[pos...] & mask) != 0
            return true
        end
        M[pos...] |= mask
        next_pos = pos .+ dir_dpos[dir + 1]
        if !checkbounds(Bool, M, next_pos...)
            return false
        end
        if M[next_pos...] != -1
            pos = next_pos
        else
            dir = (dir + 1) % 4
        end
    end
end

function count_options(file)
    M, pos = load_map(file)
    Mc = copy(M)
    c = 0
    for Idx in CartesianIndices(M)
        if M[Idx] != 0 || Tuple(Idx) == pos
            continue
        end
        Mc .= M
        Mc[Idx] = -1
        if check_loop(Mc, pos, 0)
            c += 1
        end
    end
    return c
end

@show count_options(ARGS[1])

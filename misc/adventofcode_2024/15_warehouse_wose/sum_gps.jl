#!/usr/bin/julia

function sum_gps(input, input_map)
    lines = readlines(input_map)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int}(undef, nrow, ncol)
    p = (0, 0)
    for row in 1:nrow
        line = lines[row]
        for col in 1:ncol
            c = line[col]
            if c == '#'
                M[row, col] = -1
            elseif c == 'O'
                M[row, col] = 1
            elseif c == '.'
                M[row, col] = 0
            elseif c == '@'
                M[row, col] = 0
                p = (row, col)
            else
                error("Unknown character $c")
            end
        end
    end

    for line in eachline(input)
        for c in line
            if c == '^'
                dir = (-1, 0)
            elseif c == 'v'
                dir = (1, 0)
            elseif c == '>'
                dir = (0, 1)
            elseif c == '<'
                dir = (0, -1)
            else
                error("Unknown character $c")
            end
            newp = p .+ dir
            n = 1
            while M[(p .+ dir .* n)...] == 1
                n += 1
            end
            if M[(p .+ dir .* n)...] == -1
                continue
            end
            @assert M[(p .+ dir .* n)...] == 0
            if n != 1
                M[newp...] = 0
                M[(p .+ dir .* n)...] = 1
            end
            p = newp
        end
    end

    s = 0
    for row in 1:nrow
        for col in 1:ncol
            if M[row, col] == 1
                s += (row - 1) * 100 + (col - 1)
            end
        end
    end
    return s
end

@show sum_gps(ARGS[1], ARGS[2])

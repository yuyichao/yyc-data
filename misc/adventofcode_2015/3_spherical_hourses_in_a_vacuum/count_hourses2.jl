#!/usr/bin/julia

function count_houses(file)
    pos = (0, 0)
    pos2 = (0, 0)
    visited = Set((pos,))
    for c in read(file, String)
        if c == '^'
            pos = pos .+ (0, 1)
        elseif c == 'v'
            pos = pos .+ (0, -1)
        elseif c == '<'
            pos = pos .+ (-1, 0)
        elseif c == '>'
            pos = pos .+ (1, 0)
        end
        push!(visited, pos)
        pos, pos2 = pos2, pos
    end
    return length(visited)
end

@show count_houses(ARGS[1])

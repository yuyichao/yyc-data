#!/usr/bin/julia

function find_distance(file)
    pos = (0, 0)
    dir = (0, 1)
    visited = Set((pos,))
    for inst in split(read(file, String), ',')
        m = match(r"([RL])(\d+)", inst)
        if m[1] == "R"
            dir = (dir[2], -dir[1])
        else
            dir = (-dir[2], dir[1])
        end
        dist = parse(Int, m[2])
        for i in 1:dist
            pos = pos .+ dir
            if pos in visited
                return sum(abs.(pos))
            end
            push!(visited, pos)
        end
    end
end

@show find_distance(ARGS[1])

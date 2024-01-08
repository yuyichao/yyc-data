#!/usr/bin/julia

function find_distance(file)
    pos = (0, 0)
    dir = (0, 1)
    for inst in split(read(file, String), ',')
        m = match(r"([RL])(\d+)", inst)
        if m[1] == "R"
            dir = (dir[2], -dir[1])
        else
            dir = (-dir[2], dir[1])
        end
        pos = pos .+ parse(Int, m[2]) .* dir
    end
    return sum(abs.(pos))
end

@show find_distance(ARGS[1])

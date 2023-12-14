#!/usr/bin/julia

function count_wins(file)
    lines = eachline(file)
    times = parse.(Int, split(split(first(lines), ":")[2]))
    distances = parse.(Int, split(split(first(lines), ":")[2]))

    ct = 1
    for (t, d) in zip(times, distances)
        c = 0
        for i in 1:t - 1
            if i * (t - i) > d
                c += 1
            end
        end
        ct *= c
    end
    return ct
end

@show count_wins(ARGS[1])

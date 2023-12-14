#!/usr/bin/julia

function count_wins(file)
    lines = eachline(file)
    time = parse(Int, replace(split(first(lines), ":")[2], " "=>""))
    distance = parse(Int, replace(split(first(lines), ":")[2], " "=>""))

    c = 0
    for i in 1:time - 1
        if i * (time - i) > distance
            c += 1
        end
    end
    return c
end

@show count_wins(ARGS[1])

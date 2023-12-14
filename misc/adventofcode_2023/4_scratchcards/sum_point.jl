#!/usr/bin/julia

function sum_point(file)
    s = 0
    for line in eachline(file)
        _, results = split(line, ":")
        win, mine = split(results, "|")
        win = Set(split(win))
        n = 0
        for v in split(mine)
            if v in win
                n += 1
            end
        end
        if n >= 1
            s += 2^(n - 1)
        end
    end
    return s
end

@show sum_point(ARGS[1])

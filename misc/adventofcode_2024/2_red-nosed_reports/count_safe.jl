#!/usr/bin/julia

function count_safe(file)
    c = 0
    for line in eachline(file)
        levels = parse.(Int, split(line))
        nl = length(levels)
        d = levels[2] - levels[1]
        if d > 3 || d < -3 || d == 0
            continue
        end
        safe = true
        for i in 3:nl
            d′ = levels[i] - levels[i - 1]
            if d′ > 3 || d′ < -3 || d′ == 0 || d′ * d < 0
                safe = false
                break
            end
        end
        if safe
            c += 1
        end
    end
    return c
end

@show count_safe(ARGS[1])

#!/usr/bin/julia

function is_safe(levels)
    nl = length(levels)
    d = levels[2] - levels[1]
    if d > 3 || d < -3 || d == 0
        return false
    end
    for i in 3:nl
        d′ = levels[i] - levels[i - 1]
        if d′ > 3 || d′ < -3 || d′ == 0 || d′ * d < 0
            return false
        end
    end
    return true
end

function count_safe(file)
    c = 0
    for line in eachline(file)
        levels = parse.(Int, split(line))
        nl = length(levels)
        if is_safe(levels)
            c += 1
            continue
        end
        for i in 1:nl
            l′ = [levels[j] for j in 1:nl if j != i]
            if is_safe(l′)
                c += 1
                break
            end
        end
    end
    return c
end

@show count_safe(ARGS[1])

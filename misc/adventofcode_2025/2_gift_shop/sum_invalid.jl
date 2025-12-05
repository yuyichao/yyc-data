#!/usr/bin/julia

function sum_invalid(file)
    s = 0
    for rs in split(read(file, String), ",")
        lbs, ubs = split(rs, "-")
        for v in parse(Int, lbs):parse(Int, ubs)
            vs = string(v)
            n = length(vs)
            if n % 2 == 0 && vs[1:n ÷ 2] == vs[n ÷ 2 + 1:end]
                s += v
            end
        end
    end
    return s
end

@show sum_invalid(ARGS[1])

#!/usr/bin/julia

function sum_invalid(file)
    s = 0
    for rs in split(read(file, String), ",")
        lbs, ubs = split(rs, "-")
        for v in parse(Int, lbs):parse(Int, ubs)
            if match(r"^(\d+)\1+$", string(v)) !== nothing
                s += v
            end
        end
    end
    return s
end

@show sum_invalid(ARGS[1])

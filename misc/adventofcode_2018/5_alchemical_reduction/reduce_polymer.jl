#!/usr/bin/julia

function reduce_polymer(file)
    p = collect(strip(read(file, String)))
    i = 1
    while i < length(p)
        c1 = p[i]
        c2 = p[i + 1]
        if abs(c1 - c2) == abs('A' - 'a')
            deleteat!(p, i:i + 1)
            i = max(i - 1, 1)
        else
            i += 1
        end
    end
    return length(p)
end

@show reduce_polymer(ARGS[1])

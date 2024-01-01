#!/usr/bin/julia

function reduce_polymer(p)
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

function reduce_polymer2(file)
    p = collect(strip(read(file, String)))
    return minimum(reduce_polymer(filter(c->(c != c0 && c != c0 + ('a' - 'A')), p))
                   for c0 in 'A':'Z')
end

@show reduce_polymer2(ARGS[1])

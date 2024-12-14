#!/usr/bin/julia

function calc_safe(file)
    M = zeros(Int, 101, 103)
    for line in eachline(file)
        m = match(r"p=(\d*),(\d*) v=(-?\d*),(-?\d*)", line)
        p = parse(Int, m[1]), parse(Int, m[2])
        v = parse(Int, m[3]), parse(Int, m[4])
        pn = p .+ v .* 100
        x = mod(pn[1], 101) + 1
        y = mod(pn[2], 103) + 1
        M[x, y] += 1
    end

    s = 1
    for p0 in ((1, 1), (52, 1), (1, 53), (52, 53))
        s *= sum(M[p0[1]:p0[1] + 49, p0[2]:p0[2] + 50])
    end
    return s
end

@show calc_safe(ARGS[1])

#!/usr/bin/julia

function parse_robot(file)
    r = NTuple{2,NTuple{2,Int}}[]
    for line in eachline(file)
        m = match(r"p=(\d*),(\d*) v=(-?\d*),(-?\d*)", line)
        p = parse(Int, m[1]), parse(Int, m[2])
        v = parse(Int, m[3]), parse(Int, m[4])
        push!(r, (p, v))
    end
    return r
end

const rs = parse_robot(ARGS[1])
const M = zeros(Float64, 101 * 101, 103 * 103)

function print_map(M, rs, n)
    M .= 0
    for (p, v) in rs
        pn = p .+ v .* n
        x = mod(pn[1], 101) + 1
        y = mod(pn[2], 103) + 1
        M[x, y] = 1
    end
end

for i in 1:101
    for j in 1:103
        idx = (i - 1) * 103 + j
        print_map(@view(M[(i - 1) * 101 + 1:(i - 1) * 101 + 101,
                          (j - 1) * 103 + 1:(j - 1) * 103 + 103]), rs, idx)
    end
end

println("Generated")

using Images
save("all.png", colorview(Gray, M))

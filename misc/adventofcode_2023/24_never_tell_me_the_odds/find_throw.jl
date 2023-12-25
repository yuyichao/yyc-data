#!/usr/bin/julia

using NLsolve

function solve(p1, v1, p2, v2, p3, v3, init)
    x1, y1, z1 = p1
    vx1, vy1, vz1 = v1
    x2, y2, z2 = p2
    vx2, vy2, vz2 = v2
    x3, y3, z3 = p3
    vx3, vy3, vz3 = v3
    function f(X)
        x0, y0, z0, vx0, vy0, vz0 = X
        return [(x1 - x0) * (vy1 - vy0) - (y1 - y0) * (vx1 - vx0),
                (y1 - y0) * (vz1 - vz0) - (z1 - z0) * (vy1 - vy0),
                (x2 - x0) * (vy2 - vy0) - (y2 - y0) * (vx2 - vx0),
                (y2 - y0) * (vz2 - vz0) - (z2 - z0) * (vy2 - vy0),
                (x3 - x0) * (vy3 - vy0) - (y3 - y0) * (vx3 - vx0),
                (y3 - y0) * (vz3 - vz0) - (z3 - z0) * (vy3 - vy0)]
    end
    # init = Float64[x1, y1, z1, vx1, vy1, vz1]
    @show f(init)
    sol1 = nlsolve(f, init)
    @show f(sol1.zero)
    return sol1
end

function find_throw(file)
    points = NTuple{2,NTuple{3,Int}}[]
    for line in eachline(file)
        m = match(r"([-0-9]+), ([-0-9]+), ([-0-9]+) @ ([-0-9]+), ([-0-9]+), ([-0-9]+)",
                  line)
        push!(points, ((parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])),
                       (parse(Int, m[4]), parse(Int, m[5]), parse(Int, m[6]))))
    end
    sol1 = solve(points[1][1], points[1][2],
                 points[2][1], points[2][2],
                 points[3][1], points[3][2],
                 BigFloat[points[10][1]..., points[10][2]...])
    return round.(BigInt, sol1.zero)
    # sol2 = solve(points[4][1], points[4][2],
    #              points[5][1], points[5][2],
    #              points[6][1], points[6][2],
    #              sol1.zero)
end

@show find_throw(ARGS[1])

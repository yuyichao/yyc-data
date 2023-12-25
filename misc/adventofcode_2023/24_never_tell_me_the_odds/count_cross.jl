#!/usr/bin/julia

function find_crossing(r1, v1, r2, v2)
    x1, y1 = r1
    vx1, vy1 = v1
    x2, y2 = r2
    vx2, vy2 = v2

    D = vx2 * vy1 - vx1 * vy2
    if D == 0
        # Ignore the overlap case
        return
    end
    t1 = (vx2 * (y2 - y1) + vy2 * (x1 - x2)) / D
    t2 = (vx1 * (y2 - y1) + vy1 * (x1 - x2)) / D
    if t1 < 0 || t2 < 0
        return
    end

    return x1 + vx1 * t1, y1 + vy1 * t1
end

in_range(v, rng) = rng[1] <= v <= rng[2]

function count_cross_within(points, xrng, yrng)
    npoints = length(points)
    c = 0
    for i in 2:npoints
        for j in 1:i - 1
            r1, v1 = points[i]
            r2, v2 = points[j]
            cross = find_crossing(r1, v1, r2, v2)
            if cross === nothing
                continue
            end
            xc, yc = cross
            if in_range(xc, xrng) && in_range(yc, yrng)
                c += 1
            end
        end
    end
    return c
end

function count_cross(file)
    points = NTuple{2,NTuple{2,Int}}[]
    for line in eachline(file)
        m = match(r"([-0-9]+), ([-0-9]+), ([-0-9]+) @ ([-0-9]+), ([-0-9]+), ([-0-9]+)",
                  line)
        push!(points, ((parse(Int, m[1]), parse(Int, m[2])),
                       (parse(Int, m[4]), parse(Int, m[5]))))
    end
    return count_cross_within(points, (200000000000000, 400000000000000),
                              (200000000000000, 400000000000000))
end

@show count_cross(ARGS[1])

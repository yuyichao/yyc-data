#!/usr/bin/julia

function count_overlap(file)
    cross_count = Dict{NTuple{2,Int},Int}()
    for line in eachline(file)
        m = match(r"(\d+),(\d+) -> (\d+),(\d+)", line)
        x1, y1 = parse(Int, m[1]), parse(Int, m[2])
        x2, y2 = parse(Int, m[3]), parse(Int, m[4])
        if x1 == x2
            for y in min(y1, y2):max(y1, y2)
                cross_count[(x1, y)] = get(cross_count, (x1, y), 0) + 1
            end
        elseif y1 == y2
            for x in min(x1, x2):max(x1, x2)
                cross_count[(x, y1)] = get(cross_count, (x, y1), 0) + 1
            end
        else
            @assert abs(y1 - y2) == abs(x1 - x2)
            n = abs(x1 - x2) + 1
            dx = sign(x2 - x1)
            dy = sign(y2 - y1)
            for i in 1:n
                x = x1 + dx * (i - 1)
                y = y1 + dy * (i - 1)
                cross_count[(x, y)] = get(cross_count, (x, y), 0) + 1
            end
        end
    end
    return count(>(1), values(cross_count))
end

@show count_overlap(ARGS[1])

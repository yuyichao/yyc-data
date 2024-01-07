#!/usr/bin/julia

function trace_path(line)
    x, y = 0, 0
    points = Set{NTuple{2,Int}}()
    for inst in split(line, ',')
        m = match(r"([RULD])(\d+)", inst)
        d = parse(Int, m[2])
        if m[1] == "R"
            for i in 1:d
                push!(points, (x + i, y))
            end
            x += d
        elseif m[1] == "L"
            for i in 1:d
                push!(points, (x - i, y))
            end
            x -= d
        elseif m[1] == "U"
            for i in 1:d
                push!(points, (x, y + i))
            end
            y += d
        elseif m[1] == "D"
            for i in 1:d
                push!(points, (x, y - i))
            end
            y -= d
        end
    end
    return points
end

function min_crossing(file)
    lines = eachline(file)
    path1 = trace_path(first(lines))
    path2 = trace_path(first(lines))
    crossings = intersect(path1, path2)
    minimum(abs(x) + abs(y) for (x, y) in crossings)
end

@show min_crossing(ARGS[1])

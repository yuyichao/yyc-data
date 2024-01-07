#!/usr/bin/julia

function trace_path(line)
    x = Ref(0)
    y = Ref(0)
    step = Ref(0)
    points = Dict{NTuple{2,Int},Int}()

    function take_step(dx, dy)
        step[] += 1
        x[] += dx
        y[] += dy
        get!(points, (x[], y[]), step[])
        return
    end

    for inst in split(line, ',')
        m = match(r"([RULD])(\d+)", inst)
        d = parse(Int, m[2])
        if m[1] == "R"
            for i in 1:d
                take_step(1, 0)
            end
        elseif m[1] == "L"
            for i in 1:d
                take_step(-1, 0)
            end
        elseif m[1] == "U"
            for i in 1:d
                take_step(0, 1)
            end
        elseif m[1] == "D"
            for i in 1:d
                take_step(0, -1)
            end
        end
    end
    return points
end

function min_crossing(file)
    lines = eachline(file)
    path1 = trace_path(first(lines))
    path2 = trace_path(first(lines))
    crossings = intersect(keys(path1), keys(path2))
    minimum(path1[(x, y)] + path2[(x, y)] for (x, y) in crossings)
end

@show min_crossing(ARGS[1])

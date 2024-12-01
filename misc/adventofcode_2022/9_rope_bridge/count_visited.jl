#!/usr/bin/julia

function move_head(hpos, tpos, dir, s)
    if dir == "U"
        hpos = hpos .+ (0, 1)
    elseif dir == "D"
        hpos = hpos .+ (0, -1)
    elseif dir == "L"
        hpos = hpos .+ (-1, 0)
    elseif dir == "R"
        hpos = hpos .+ (1, 0)
    else
        error("Invalid direction $dir")
    end
    diff = tpos .- hpos
    if diff[1] > 1
        diff = (1, 0)
    elseif diff[1] < -1
        diff = (-1, 0)
    elseif diff[2] > 1
        diff = (0, 1)
    elseif diff[2] < -1
        diff = (0, -1)
    end
    tpos = hpos .+ diff
    push!(s, tpos)
    return hpos, tpos
end

function count_visited(file)
    tpos = (0, 0)
    hpos = (0, 0)
    s = Set((tpos,))
    for line in eachline(file)
        dir, _step = split(line)
        for _ in 1:parse(Int, _step)
            hpos, tpos = move_head(hpos, tpos, dir, s)
        end
    end
    return length(s)
end

@show count_visited(ARGS[1])

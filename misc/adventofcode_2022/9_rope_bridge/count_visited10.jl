#!/usr/bin/julia

function move_head!(rope, dir)
    hpos = rope[1]
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
    rope[1] = hpos
    return
end

function move_node!(rope, idx)
    prev = rope[idx - 1]
    cur = rope[idx]
    diff = cur .- prev
    if abs(diff[1]) <= 1 && abs(diff[2]) <= 1
        return
    end
    diff = diff .- sign.(diff)
    rope[idx] = prev .+ diff
    return
end

function move_all_node!(rope)
    for i in 2:length(rope)
        move_node!(rope, i)
    end
end

function count_visited(file)
    rope = [(0, 0) for _ in 1:10]
    s = Set((rope[end],))
    for line in eachline(file)
        dir, _step = split(line)
        for _ in 1:parse(Int, _step)
            move_head!(rope, dir)
            move_all_node!(rope)
            push!(s, rope[end])
        end
    end
    return length(s)
end

@show count_visited(ARGS[1])

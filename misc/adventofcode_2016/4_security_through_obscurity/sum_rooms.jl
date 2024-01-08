#!/usr/bin/julia

using DataStructures

function real_room_id(line)
    m = match(r"(.*)-(\d+)\[(.....)\]", line)
    a = counter(m[1])
    reset!(a, '-')
    counts = collect(a)
    sort!(counts, by=x->(-x[2], x[1]))
    if String([counts[i][1] for i in 1:5]) == m[3]
        return parse(Int, m[2])
    end
    return 0
end

@show sum(real_room_id, eachline(ARGS[1]))

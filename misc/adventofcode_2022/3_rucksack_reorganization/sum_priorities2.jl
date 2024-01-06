#!/usr/bin/julia

function get_priority(line1, line2, line3)
    p = intersect(line1, line2, line3)
    @assert length(p) == 1
    p = p[1]
    if 'a' <= p <= 'z'
        return p - 'a' + 1
    else
        @assert 'A' <= p <= 'Z'
        return p - 'A' + 27
    end
end

function sum_priorities(file)
    lines = readlines(file)
    s = 0
    for i in 1:3:length(lines)
        s += get_priority(lines[i], lines[i + 1], lines[i + 2])
    end
    return s
end

@show sum_priorities(ARGS[1])

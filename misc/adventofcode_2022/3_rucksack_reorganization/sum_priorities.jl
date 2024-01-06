#!/usr/bin/julia

function get_priority(line)
    len = length(line)
    p = intersect(line[1:end ÷ 2], line[end ÷ 2 + 1:end])
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
    return sum(get_priority(line) for line in eachline(file))
end

@show sum_priorities(ARGS[1])

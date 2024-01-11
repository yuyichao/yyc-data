#!/usr/bin/julia

function length_diff(file)
    s = 0
    for line in eachline(file)
        s += length(repr(line)) - length(line)
    end
    return s
end

@show length_diff(ARGS[1])

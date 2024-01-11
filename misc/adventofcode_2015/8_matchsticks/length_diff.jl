#!/usr/bin/julia

function length_diff(file)
    s = 0
    for line in eachline(file)
        str = Meta.parse(line)::String
        s += length(line) - length(str)
    end
    return s
end

@show length_diff(ARGS[1])

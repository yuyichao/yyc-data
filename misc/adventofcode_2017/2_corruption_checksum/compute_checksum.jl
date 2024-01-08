#!/usr/bin/julia

function compute_checksum(file)
    return -sum(-(extrema(parse.(Int, split(line)))...) for line in eachline(file))
end

@show compute_checksum(ARGS[1])

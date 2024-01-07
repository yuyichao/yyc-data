#!/usr/bin/julia

function calc_frequency(file)
    sum(parse(Int, line) for line in eachline(file))
end

@show calc_frequency(ARGS[1])

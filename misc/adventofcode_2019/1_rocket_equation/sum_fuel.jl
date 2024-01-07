#!/usr/bin/julia

function sum_fuel(file)
    sum(parse(Int, line) ÷ 3 - 2 for line in eachline(file))
end

@show sum_fuel(ARGS[1])

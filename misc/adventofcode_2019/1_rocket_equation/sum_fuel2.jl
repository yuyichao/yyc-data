#!/usr/bin/julia

function calc_fuel(m)
    s = 0
    while true
        m = m ÷ 3 - 2
        if m <= 0
            return s
        end
        s += m
    end
end

function sum_fuel(file)
    sum(calc_fuel(parse(Int, line)) for line in eachline(file))
end

@show sum_fuel(ARGS[1])

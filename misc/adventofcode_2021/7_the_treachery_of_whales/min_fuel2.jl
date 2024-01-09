#!/usr/bin/julia

function calc_fuel(dn)
    dn = abs(dn)
    return dn * (dn + 1) ÷ 2
end

calc_all_fuel(pos, c) = sum(calc_fuel(p - c) for p in pos)

function min_fuel(file)
    pos = sort!(parse.(Int, split(read(file, String), ',')))
    return findmin(c->calc_all_fuel(pos, c), pos[1]:pos[end])
end

@show min_fuel(ARGS[1])

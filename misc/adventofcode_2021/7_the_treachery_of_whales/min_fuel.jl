#!/usr/bin/julia

function min_fuel(file)
    pos = sort!(parse.(Int, split(read(file, String), ',')))
    mid = pos[(end + 1) ÷ 2]
    return sum(abs.(pos .- mid))
end

@show min_fuel(ARGS[1])

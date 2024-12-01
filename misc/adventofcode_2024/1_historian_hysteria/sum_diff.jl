#!/usr/bin/julia

function sum_hash(file)
    p1 = Int[]
    p2 = Int[]
    for line in eachline(file)
        v1, v2 = parse.(Int, split(line))
        push!(p1, v1)
        push!(p2, v2)
    end
    sort!(p1)
    sort!(p2)
    return sum(abs.(p1 .- p2))
end

@show sum_hash(ARGS[1])

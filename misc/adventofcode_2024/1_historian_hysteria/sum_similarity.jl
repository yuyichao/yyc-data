#!/usr/bin/julia

using DataStructures

function sum_similarity(file)
    p1 = Int[]
    p2 = Accumulator{Int,Int}()
    for line in eachline(file)
        v1, v2 = parse.(Int, split(line))
        push!(p1, v1)
        push!(p2, v2)
    end
    return sum(v1 * p2[v1] for v1 in p1)
end

@show sum_similarity(ARGS[1])

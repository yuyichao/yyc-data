#!/usr/bin/julia

function count(file)
    cs = Int[]
    c = 0
    for line in eachline(file)
        if isempty(line)
            push!(cs, c)
            c = 0
        else
            c += parse(Int, line)
        end
    end
    push!(cs, c)
    sort!(cs)
    return cs[end] + cs[end - 1] + cs[end - 2]
end

@show count(ARGS[1])

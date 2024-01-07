#!/usr/bin/julia

function find_frequency(file)
    dfs = [parse(Int, line) for line in eachline(file)]
    f = 0
    fs = Set{Int}()
    for df in Iterators.cycle(dfs)
        if f in fs
            return f
        end
        push!(fs, f)
        f += df
    end
end

@show find_frequency(ARGS[1])

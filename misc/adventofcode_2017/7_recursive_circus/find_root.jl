#!/usr/bin/julia

function find_root(file)
    parents = Set{String}()
    children = Set{String}()

    for line in eachline(file)
        line = split(line, "->")
        if length(line) != 2
            continue
        end
        m = match(r"([^ ]+) \(\d+\)", line[1])
        push!(parents, m[1])
        for c in split(line[2], ',')
            push!(children, strip(c))
        end
    end

    return setdiff(parents, children)
end

@show find_root(ARGS[1])

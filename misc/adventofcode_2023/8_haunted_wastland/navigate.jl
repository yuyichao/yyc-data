#!/usr/bin/julia

function load_file(file)
    lines = eachline(file)
    direction = first(lines)
    nodes = Dict{String,NTuple{2,String}}()

    for line in lines
        if isempty(line)
            continue
        end
        m = match(r"(...) = \((...), (...)\)", line)
        nodes[m[1]] = (m[2], m[3])
    end
    return direction, nodes
end

function navigate(direction, nodes)
    s = 0
    pos = "AAA"
    for dir in Iterators.cycle(direction)
        pos = nodes[pos][dir == 'L' ? 1 : 2]
        s += 1
        if pos == "ZZZ"
            break
        end
    end
    return s
end

@show navigate(load_file(ARGS[1])...)

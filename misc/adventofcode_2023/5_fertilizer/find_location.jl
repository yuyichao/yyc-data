#!/usr/bin/julia

function get_map(map::Vector{NTuple{3,Int}}, v)
    for (n1, n2, n3) in map
        if n2 <= v <= n2 + n3 - 1
            return v - n2 + n1
        end
    end
    return v
end

function read_map(it)
    map = NTuple{3,Int}[]
    for line in it
        if isempty(line)
            break
        end
        n1, n2, n3 = parse.(Int, split(line))
        push!(map, (n1, n2, n3))
    end
    return map
end

function find_location(file)
    lines = eachline(file)

    items = Set(parse.(Int, split(split(first(lines), ":")[2])))
    name = "seed"

    while true
        line = first(lines)
        if isempty(line)
            continue
        end
        m = match(r"([^-]*)-to-([^-]*) map:", line)
        @assert m !== nothing
        @assert m[1] == name
        name = String(m[2])
        map = read_map(lines)
        items = Set(get_map(map, item) for item in items)
        if name == "location"
            break
        end
    end
    return minimum(items)
end

@show find_location(ARGS[1])

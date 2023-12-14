#!/usr/bin/julia

function convert_range(out, map::Vector{NTuple{3,Int}}, range)
    # Assuming map is sorted
    lo, hi = range
    for (n1, n2, n3) in map
        if n2 > hi
            push!(out, (lo, hi))
            return
        elseif n2 + n3 - 1 < lo
            continue
        elseif n2 + n3 - 1 >= hi
            if n2 <= lo
                push!(out, (lo, hi) .+ (n1 - n2))
            else
                push!(out, (lo, n2 - 1))
                push!(out, (n1, hi + (n1 - n2)))
            end
            return
        else
            if n2 <= lo
                push!(out, (lo + (n1 - n2), n1 + n3 - 1))
            else
                push!(out, (lo, n2 - 1))
                push!(out, (n1, n1 + n3 - 1))
            end
            lo = n2 + n3
        end
    end
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
    sort!(map, by=((n1, n2, n3),)->(n2, n3, n1))
    return map
end

function find_location(file)
    lines = eachline(file)

    items = parse.(Int, split(split(first(lines), ":")[2]))
    item_range = [(items[i], items[i] + items[i + 1] - 1) for i in 1:2:length(items)]
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
        item_range2 = NTuple{2,Int}[]
        for rng in item_range
            convert_range(item_range2, map, rng)
        end
        item_range = item_range2
        if name == "location"
            break
        end
    end
    sort!(item_range)
    return minimum(item_range)
end

@show find_location(ARGS[1])

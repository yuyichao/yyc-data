#!/usr/bin/julia

function collect_corners(instructions)
    corners = Dict{Int,Set{Int}}()
    last_dir = instructions[end][1]
    x = 0
    y = 0
    for (dir, steps) in instructions
        @assert last_dir != dir
        push!(get!(Set{Int}, corners, y), x)
        last_dir = dir
        if dir == 'R' || dir == '0'
            x += steps
        elseif dir == 'L' || dir == '2'
            x -= steps
        elseif dir == 'U' || dir == '3'
            y -= steps
        elseif dir == 'D' || dir == '1'
            y += steps
        else
            error()
        end
    end
    return corners
end

function process_line!(out_edges, in_edges, corners)
    empty!(out_edges)
    for c in in_edges
        push!(out_edges, c)
    end
    for c in corners
        if c in out_edges
            delete!(out_edges, c)
        else
            push!(out_edges, c)
        end
    end
end

function count_line_area(edges)
    edges = sort!(collect(edges))
    @assert length(edges) % 2 == 0
    a = 0
    for i in 1:2:length(edges)
        a += edges[i + 1] - edges[i] + 1
    end
    return a
end

function count_line_area(edges1, edges2)
    edges1 = sort!(collect(edges1))
    n1 = length(edges1)
    @assert n1 % 2 == 0
    edges2 = sort!(collect(edges2))
    n2 = length(edges2)
    @assert n2 % 2 == 0
    a = 0

    i1 = 1
    i2 = 1
    last_end = typemin(Int)
    while i1 <= n1 || i2 <= n2
        x11 = i1 <= n1 ? edges1[i1] : typemax(Int)
        x21 = i2 <= n2 ? edges2[i2] : typemax(Int)
        if x11 < x21
            x1 = x11
            x2 = edges1[i1 + 1]
            i1 += 2
        else
            x1 = x21
            x2 = edges2[i2 + 1]
            i2 += 2
        end
        if last_end < x2
            a += x2 - max(x1 - 1, last_end)
            last_end = x2
        end
    end
    return a
end

function count_area(file)
    instructions = Tuple{Char,Int}[]
    for line in eachline(file)
        m = match(r"^([RLUD]) ([0-9]+) \(#([0-9a-z]{5})([0-9a-z])\)", line)
        push!(instructions, (m[4][1], parse(Int, m[3], base=16)))
    end
    corners = collect_corners(instructions)
    lines = sort!(collect(keys(corners)))

    edges1 = Set{Int}()
    edges2 = Set{Int}()

    a = 0
    active_area = 0
    last_line = 0
    for l in lines
        a += active_area * (l - last_line - 1)
        corner = corners[l]
        process_line!(edges2, edges1, corner)
        a += count_line_area(edges1, edges2)
        last_line = l
        active_area = count_line_area(edges2)
        edges2, edges1 = edges1, edges2
    end
    @assert active_area == 0
    return a
end

@show count_area(ARGS[1])

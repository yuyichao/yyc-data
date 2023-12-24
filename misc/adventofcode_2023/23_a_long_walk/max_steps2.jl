#!/usr/bin/julia

const LEFT = 1
const RIGHT = 2
const UP = 3
const DOWN = 4

# function check_pos(M, visited, x, y)
#     if !checkbounds(Bool, M, y, x) || M[y, x] == -1
#         return false
#     end
#     if (x, y) in visited
#         return false
#     end
#     return true
# end

# const active_count = Ref(0)

function trace_max_steps(nodes, pos, step, end_pos, visited, max_steps)
    while true
        if max_steps[pos] < step
            max_steps[pos] = step
        end
        if end_pos == pos
            return
        end
        push!(visited, pos)

        nodeinfo = nodes[pos]

        for i in 2:length(nodeinfo.next)
            x, y, len = nodeinfo.next[i]
            if (x, y) in visited
                continue
            end
            trace_max_steps(nodes, (x, y), step + len, end_pos,
                            copy(visited), max_steps)
        end
        x, y, len = nodeinfo.next[1]
        if (x, y) in visited
            return
        end
        step += len
        pos = (x, y)
    end
end

function trace_max_steps(nodes, start_pos, end_pos)
    max_steps = Dict(pos=>0 for pos in keys(nodes))
    visited = Set{NTuple{2,Int}}()
    trace_max_steps(nodes, start_pos, 0, end_pos, visited, max_steps)
    return max_steps[end_pos]
end

struct NodeInfo
    next::Vector{NTuple{3,Int}} # (x, y, steps)
    NodeInfo() = new(NTuple{3,Int}[])
end

function find_nodes(M)
    nodes = Dict{NTuple{2,Int},NodeInfo}()
    nrows, ncols = size(M)
    for x in 1:ncols
        for y in 1:nrows
            if M[y, x] != 0
                continue
            end
            nneighbor = 0
            for (x2, y2) in ((x - 1, y), (x + 1, y),
                             (x, y - 1), (x, y + 1))
                if checkbounds(Bool, M, y2, x2) && M[y2, x2] == 0
                    nneighbor += 1
                end
            end
            if nneighbor != 0 && nneighbor != 2
                nodes[(x, y)] = NodeInfo()
            end
        end
    end
    return nodes
end

function find_next_node(M, x0, y0, x1, y1, nodes)
    steps = 1
    while true
        if haskey(nodes, (x1, y1))
            return (x1, y1, steps)
        end
        steps += 1
        found = false
        for (x2, y2) in ((x1 - 1, y1), (x1 + 1, y1),
                         (x1, y1 - 1), (x1, y1 + 1))
            if (x2, y2) == (x0, y0)
                continue
            end
            if checkbounds(Bool, M, y2, x2) && M[y2, x2] == 0
                x0, y0 = x1, y1
                x1, y1 = x2, y2
                found = true
                break
            end
        end
        @assert found
    end
end

function fill_trace_nodeinfo!(M, x, y, nodeinfo, nodes)
    for (x2, y2) in ((x - 1, y), (x + 1, y),
                     (x, y - 1), (x, y + 1))
        if checkbounds(Bool, M, y2, x2) && M[y2, x2] == 0
            push!(nodeinfo.next, find_next_node(M, x, y, x2, y2, nodes))
        end
    end
end

function max_steps(file)
    frontier = NTuple{2,Int}[]

    lines = readlines(file)
    nrows = length(lines)
    ncols = length(lines[1])

    M = Matrix{Int8}(undef, nrows, ncols)
    start_pos = (0, 0)
    end_pos = (0, 0)
    for i in 1:nrows
        line = lines[i]
        for j in 1:ncols
            c = line[j]
            if c == '#'
                M[i, j] = -1
            else
                M[i, j] = 0
                if i == 1
                    @assert start_pos == (0, 0)
                    start_pos = (j, i)
                elseif i == nrows
                    @assert end_pos == (0, 0)
                    end_pos = (j, i)
                end
            end
        end
    end

    nodes = find_nodes(M)
    @assert haskey(nodes, start_pos)
    @assert haskey(nodes, end_pos)

    for ((x, y), nodeinfo) in nodes
        fill_trace_nodeinfo!(M, x, y, nodeinfo, nodes)
    end
    trace_max_steps(nodes, start_pos, end_pos)
end

@show max_steps(ARGS[1])

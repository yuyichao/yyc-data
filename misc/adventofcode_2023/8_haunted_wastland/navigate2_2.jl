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
    direction = [dir == 'R' for dir in direction]
    node_names = collect(keys(nodes))
    node_ids = Dict{String,Int}()
    for (i, name) in enumerate(node_names)
        node_ids[name] = i
    end
    left_next = [node_ids[nodes[name][1]] for name in node_names]
    right_next = [node_ids[nodes[name][2]] for name in node_names]
    is_end = [name[3] == 'Z' for name in node_names]
    starts = [node_ids[name] for name in node_names if name[3] == 'A']

    return direction, left_next, right_next, is_end, starts
end

struct NodeInfo
    next::Int
    end_pos::Set{Int}
end

struct RoundInfo
    round_len::Int
    nodes::Vector{NodeInfo}
end

function to_round(direction, left_next, right_next, is_end)
    nnodes = length(left_next)
    round_len = length(direction)
    round_next = [1:nnodes;]
    end_pos = [Set{Int}() for i in 1:nnodes]
    @inbounds for i in 1:round_len
        dir = direction[i]
        dir_next = dir ? right_next : left_next
        for j in 1:nnodes
            pos = round_next[j]
            if is_end[pos]
                push!(end_pos[j], i - 1)
            end
            round_next[j] = dir_next[pos]
        end
    end
    return RoundInfo(round_len, [NodeInfo(round_next[i], end_pos[i]) for i in 1:nnodes])
end

function navigate(file)
    direction, left_next, right_next, is_end, starts = load_file(file)
    round = to_round(direction, left_next, right_next, is_end)

    c = 1

    for s0 in starts
        seen = Dict{Int,Int}()

        end_pos = Int[]
        p = s0
        r = 0
        while true
            if haskey(seen, p)
                break
            end
            seen[p] = r
            info = round.nodes[p]
            if 0 in info.end_pos
                push!(end_pos, r)
            end
            r += 1
            p = info.next
        end
        @assert seen[p] == 1
        @assert length(end_pos) == 1
        @assert end_pos[1] + 1 == r
        c = lcm(c, end_pos[1])
    end
    return c * round.round_len
end

@show navigate(ARGS[1])

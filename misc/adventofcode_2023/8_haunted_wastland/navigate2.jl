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

function double_up(round::RoundInfo, max_end=Ref(0), min_end=Ref(typemax(Int)))
    return RoundInfo(round.round_len * 2,
                     [begin
                          next_info = round.nodes[info.next]
                          end_pos = copy(info.end_pos)
                          for pos in next_info.end_pos
                              push!(end_pos, pos + round.round_len)
                          end
                          max_end[] = max(max_end[], length(end_pos))
                          min_end[] = min(min_end[], length(end_pos))
                          NodeInfo(next_info.next, end_pos)
                      end for info in round.nodes])
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

function navigate_round(round::RoundInfo, starts)
    npos = length(starts)
    poss = copy(starts)
    end_pos = Vector{Set{Int}}(undef, npos)
    nrounds = 0
    real_end = Set{Int}()
    @inbounds while true
        min_len = typemax(Int)
        min_idx = 0
        for i in 1:npos
            pos = poss[i]
            node_info = round.nodes[pos]
            nends = length(node_info.end_pos)
            end_pos[i] = node_info.end_pos
            poss[i] = node_info.next
            if nends <= min_len
                min_len = nends
                min_idx = i
            end
        end
        for p in end_pos[min_idx]
            all_found = true
            for i in 1:npos
                if i == min_idx
                    continue
                end
                if !(p in end_pos[i])
                    all_found = false
                    break
                end
            end
            if all_found
                push!(real_end, p)
            end
        end
        if !isempty(real_end)
            return nrounds * round.round_len + minimum(real_end)
        end
        nrounds += 1
    end
end

function navigate(file)
    direction, left_next, right_next, is_end, starts = load_file(file)
    round = to_round(direction, left_next, right_next, is_end)
    while true
        max_end_count = Ref(0)
        min_end_count = Ref(typemax(Int))
        round = double_up(round, max_end_count, min_end_count)
        @show min_end_count[], max_end_count[], round.round_len
        if max_end_count[] > 65536
            break
        end
    end
    return navigate_round(round, starts)
end

@show navigate(ARGS[1])

#!/usr/bin/julia

function count_party(file)
    conn = Dict{String,Set{String}}()
    nodes = Set{String}()
    for line in eachline(file)
        node1, node2 = split(line, "-")
        push!(nodes, node1)
        push!(nodes, node2)
        push!(get!(()->Set{String}(), conn, node1), node2)
        push!(get!(()->Set{String}(), conn, node2), node1)
    end
    parties = Set{Set{String}}()
    for (node1, nodes) in conn
        for node2 in nodes
            conn2 = conn[node2]
            for node3 in nodes
                if node3 == node2
                    continue
                end
                if node3 in conn2
                    push!(parties, Set((node1, node2, node3)))
                end
            end
        end
    end
    s = 0
    for nodes in parties
        for node in nodes
            if node[1] == 't'
                s += 1
                break
            end
        end
    end
    return s
end

@show count_party(ARGS[1])

#!/usr/bin/julia

function check_party(party, conns)
    if isempty(party)
        return true
    end
    for node in party
        conn = conns[node]
        for node2 in party
            if node2 == node
                continue
            end
            if !(node2 in conn)
                return false
            end
        end
    end
    return true
end

function delete_node!(conns, node)
    delete!(conns, node)
    for (node2, conn2) in conns
        delete!(conn2, node)
    end
end

function find_max_node(conns, node1, max_party_ref)
    conn1 = collect(conns[node1])
    nconn1 = length(conn1)
    if nconn1 + 1 < length(max_party_ref[])
        return
    end

    for i in 0:(2^nconn1 - 1)
        party = Set((node1,))
        for j in 1:nconn1
            if (i >> (j - 1)) & 1 != 0
                push!(party, conn1[j])
            end
        end
        if length(party) < length(max_party_ref[])
            continue
        end
        if check_party(party, conns)
            @show length(party), party
            max_party_ref[] = party
        end
    end
end

function max_party(file)
    conns = Dict{String,Set{String}}()
    for line in eachline(file)
        node1, node2 = split(line, "-")
        push!(get!(()->Set{String}(), conns, node1), node2)
        push!(get!(()->Set{String}(), conns, node2), node1)
    end

    max_party_ref = Ref(Set{String}())
    while !isempty(conns)
        node, _ = first(conns)
        find_max_node(conns, node, max_party_ref)
        delete_node!(conns, node)
    end
    return join(sort(collect(max_party_ref[])), ",")
end

@show max_party(ARGS[1])

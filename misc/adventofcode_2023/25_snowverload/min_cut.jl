#!/usr/bin/julia

struct Node
    nodes::Set{String}
    edges::Dict{Node,Int}
end

function load_graph(file)
    nodes = Dict{String,Node}()
    function find_node(name)
        return get!(()->Node(Set{String}((name,)), Dict{Node,Int}()), nodes, name)
    end
    for line in eachline(file)
        m = match(r"([a-z]+): (.*)", line)
        name = m[1]
        edges = split(m[2])
        node = find_node(name)
        for edge in edges
            node2 = find_node(edge)
            node.edges[node2] = 1
            node2.edges[node] = 1
        end
    end
    return Set(values(nodes))
end

function merge_node(nodes, node1, node2)
    delete!(nodes, node1)
    delete!(nodes, node2)
    node = Node(union(node1.nodes, node2.nodes), Dict{Node,Int}())
    push!(nodes, node)
    for (edge, w) in node1.edges
        if edge === node2
            continue
        end
        new_w = get(node.edges, edge, 0) + w
        node.edges[edge] = new_w
        delete!(edge.edges, node1)
        edge.edges[node] = new_w
    end
    for (edge, w) in node2.edges
        if edge === node1
            continue
        end
        new_w = get(node.edges, edge, 0) + w
        node.edges[edge] = new_w
        delete!(edge.edges, node2)
        edge.edges[node] = new_w
    end
    return node
end

function compute_weight(nodes, node::Node)
    sw = 0
    for (edge, w) in node.edges
        if edge in nodes
            sw += w
        end
    end
    return sw
end

function compute_weight(node::Node)
    sw = 0
    for (edge, w) in node.edges
        sw += w
    end
    return sw
end

function count_node(node::Node)
    return length(node.nodes)
end

function count_node(nodes)
    return sum(count_node(node) for node in nodes)
end

function find_min_st_cut(nodes)
    @assert length(nodes) > 1

    if length(nodes) == 2
        node1, node2 = nodes
        c1 = count_node(node1)
        c2 = count_node(node2)
        cut_weight = compute_weight(node1)
        merge_node(nodes, node1, node2)

        return c1, c2, cut_weight, nodes
    end

    a = Set((pop!(nodes),))

    while true
        max_node = first(nodes)
        max_weight = compute_weight(a, max_node)
        for node in nodes
            if node === max_node
                continue
            end
            weight = compute_weight(a, node)
            if weight > max_weight
                max_node = node
                max_weight = weight
            end
        end
        delete!(nodes, max_node)
        push!(a, max_node)
        if length(nodes) > 1
            continue
        end
        t = max_node
        s = first(nodes)
        c1 = count_node(a)
        c2 = count_node(s)
        cut_weight = compute_weight(s)
        push!(a, s)
        new_node = merge_node(a, s, t)
        @show new_node.nodes
        return c1, c2, cut_weight, a
    end
end

function find_min_cut(nodes)
    min_c1 = 0
    min_c2 = 0
    min_cut = typemax(Int)
    while length(nodes) > 1
        c1, c2, cut, nodes = find_min_st_cut(nodes)
        @show c1, c2, cut, length(nodes)
        if cut < min_cut
            min_c1 = c1
            min_c2 = c2
            min_cut = cut
        end
    end
    return min_c1, min_c2, min_cut
end

function min_cut(file)
    nodes = load_graph(file)
    return find_min_cut(nodes)
end

@show min_cut(ARGS[1])

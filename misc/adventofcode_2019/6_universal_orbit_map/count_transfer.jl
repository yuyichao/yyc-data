#!/usr/bin/julia

mutable struct Object
    const name::String
    const direct::Vector{Object}
    depth::Int
    parent::Object
    Object(name) = new(name, Object[], 0)
end

get_object(objs, name) = get!(()->Object(name), objs, name)

function mark_depth(obj, depth)
    obj.depth = depth
    for o in obj.direct
        mark_depth(o, depth + 1)
    end
end

function find_merge_base(obj1, obj2)
    objs = Set{Object}()
    obj = obj1
    while true
        push!(objs, obj)
        if !isdefined(obj, :parent)
            break
        end
        obj = obj.parent
    end
    while true
        if obj2 in objs
            return obj2
        end
        obj2 = obj2.parent
    end
end

function count_transfer(file)
    objs = Dict{String,Object}()
    for line in eachline(file)
        c, o = split(line, ')')
        c = get_object(objs, c)
        o = get_object(objs, o)
        o.parent = c
        push!(c.direct, o)
    end
    mark_depth(objs["COM"], 0)
    you = objs["YOU"]
    san = objs["SAN"]
    base = find_merge_base(you, san)
    return (you.depth + san.depth - base.depth * 2 - 2)
end

@show count_transfer(ARGS[1])

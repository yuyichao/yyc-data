#!/usr/bin/julia

mutable struct Object
    const name::String
    const direct::Vector{Object}
    depth::Int
    Object(name) = new(name, Object[], 0)
end

get_object(objs, name) = get!(()->Object(name), objs, name)

function mark_depth(obj, depth)
    obj.depth = depth
    for o in obj.direct
        mark_depth(o, depth + 1)
    end
end

function count_orbit(file)
    objs = Dict{String,Object}()
    for line in eachline(file)
        c, o = split(line, ')')
        c = get_object(objs, c)
        o = get_object(objs, o)
        push!(c.direct, o)
    end
    mark_depth(objs["COM"], 0)
    s = 0
    for (name, obj) in objs
        @assert obj.name == "COM" || obj.depth > 0
        s += obj.depth
    end
    return s
end

@show count_orbit(ARGS[1])

#!/usr/bin/julia

using EzXML

svg = readxml(ARGS[1])
svgroot = root(svg)
g10 = firstelement(svgroot)

function fix_xy(s)
    xy = split(s, ",")
    @assert length(xy) == 2
    xy[2] = string(-parse(Float64, xy[2]))
    return join(xy, ",")
end

function fix_y(s)
    return string(-parse(Float64, s))
end

function fix_path(d)
    res = String[]
    items = split(d, " ")
    i = 1
    nitems = length(items)
    while i <= nitems
        cmd = items[i]
        push!(res, cmd)
        i += 1
        if cmd == "z"
            @assert i == nitems + 1
            break
        elseif cmd == "m"
            while !isletter(items[i][1])
                push!(res, fix_xy(items[i]))
                i += 1
            end
        elseif cmd == "M"
            while !isletter(items[i][1])
                push!(res, fix_xy(items[i]))
                i += 1
            end
        elseif cmd == "h"
            push!(res, items[i])
            i += 1
        elseif cmd == "H"
            push!(res, items[i])
            i += 1
        elseif cmd == "v"
            push!(res, fix_y(items[i]))
            i += 1
        elseif cmd == "V"
            push!(res, fix_y(items[i]))
            i += 1
        elseif cmd == "l"
            while !isletter(items[i][1])
                push!(res, fix_xy(items[i]))
                i += 1
            end
        elseif cmd == "L"
            while !isletter(items[i][1])
                push!(res, fix_xy(items[i]))
                i += 1
            end
        else
            error("Unknown command $cmd in $items")
        end
    end
    return join(res, " ")
end

for g in elements(g10)
    g.name == "g" || continue
    haskey(g, "transform") || continue
    @assert g["transform"] == "scale(1,-1)"
    id = g["id"]
    paths = elements(g)
    @assert length(paths) == 1
    path = paths[1]
    path["id"] = id
    unlink!(path)
    linkprev!(g, path)
    unlink!(g)
    path["d"] = fix_path(path["d"])
end
println(svgroot)

#!/usr/bin/julia

struct Color
    name::String
    contained_in::Set{Color}
    contain::Dict{Color,Int}
    Color(name) = new(name, Set{Color}(), Dict{Color,Int}())
end

function load_colors(file)
    colors = Dict{String,Color}()
    get_color(name) = get!(()->Color(name), colors, name)

    for line in eachline(file)
        subject, objects = split(line, " contain ", limit=2)
        ms = match(r"(.+) +bag", subject)
        subject = get_color(strip(ms[1]))
        if objects == "no other bags."
            continue
        end
        objects = split(strip(objects, ('.',)), ", ")
        for object in objects
            mo = match(r"(\d+) (.+) +bag", object)
            object = get_color(strip(mo[2]))
            subject.contain[object] = parse(Int, mo[1])
            push!(object.contained_in, subject)
        end
    end

    return colors
end

function scan_contains(color, res)
    for c in color.contained_in
        push!(res, c)
        scan_contains(c, res)
    end
    return res
end

function count_contains(file)
    colors = load_colors(file)
    return length(scan_contains(colors["shiny gold"], Set{Color}()))
end

@show count_contains(ARGS[1])

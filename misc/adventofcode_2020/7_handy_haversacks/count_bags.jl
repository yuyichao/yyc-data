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

function get_bag_count(res, color)
    return get!(res, color) do
        if isempty(color.contain)
            return 0
        end
        return sum((get_bag_count(res, c) + 1) * n for (c, n) in color.contain)
    end
end

function count_bags(file)
    colors = load_colors(file)
    return get_bag_count(Dict{Color,Int}(), colors["shiny gold"])
end

@show count_bags(ARGS[1])

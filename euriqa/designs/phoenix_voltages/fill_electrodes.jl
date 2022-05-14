#!/usr/bin/julia

using PyPlot
using Printf

const template = read(joinpath(@__DIR__, "template.svg.in"), String)

# Very dumb template engine
function fill_template(context, def)
    replace(template, r"\${\w*}"=>s->get(context, s[3:end - 1], def))
end

function val_to_color(val)
    @sprintf "%02x" clamp(round(Int, val * 0xff), UInt8)
end

function default_fill_map(v, cmap=get_cmap("RdBu_r"))
    rgba = cmap(v / 20 + 0.5)
    rgb = rgba[1:3] .* rgba[4]
    return "#" * val_to_color(rgb[1]) * val_to_color(rgb[2]) * val_to_color(rgb[3])
end

function fill_electrode(vals, fill_map=default_fill_map, def_fill="none")
    context = Dict{String,String}()
    for (key, value) in vals
        context["fill_$key"] = fill_map(value)
    end
    return fill_template(context, def_fill)
end

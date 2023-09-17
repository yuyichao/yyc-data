#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("read_data.jl")

function count_cycles(bin)
    data = bin.data
    minv, maxv = extrema(data)
    top_thresh = 0.8 * maxv + 0.2 * minv
    bottom_thresh = 0.2 * maxv + 0.8 * minv

    # 1: raise through top_thresh
    # 2: lower through bottom_thresh
    state = data[1] >= top_thresh ? 2 : 1
    raise_count = 0
    for v in data
        if state == 1
            if v > top_thresh
                raise_count += 1
                state = 2
            end
        else
            if v < bottom_thresh
                state = 1
            end
        end
    end
    spacing = length(data) / raise_count
    period = spacing * bin.dx
    return raise_count, period, Int(maxv) - Int(minv)
end

for f in readdir(joinpath(@__DIR__, "data"), join=true)
    bin = read_bin_compressed(f)
    @show f, count_cycles(bin)
end

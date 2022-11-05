#!/usr/bin/julia

using HDF5

load_data(fname) = h5open(fname) do io
    return (time=read(io, "time"), value=read(io, "value"),
            tscale=read(io, "tscale"), vscale=read(io, "vscale"))
end

function shrink_block(block, r_start, r_end)
    i = first(block)
    l = last(block)
    new_i = round(Int, i * (1 - r_start) + l * r_start)
    new_l = round(Int, i * (1 - r_end) + l * r_end)
    return new_i:new_l
end

function find_blocks(data)
    threshold = 20
    min_size = 50
    end_threshold = 100

    blocks = Tuple{Int,Int}[]
    in_block = false
    block_start = 0
    last_above = 0
    for (i, v) in enumerate(data.value)
        if v < threshold
            if !in_block
                continue
            end
            if i > last_above + end_threshold
                if last_above - block_start > min_size
                    push!(blocks, shrink_block(block_start:last_above, 0.1, 0.9))
                end
                in_block = false
            end
            continue
        end
        last_above = i
        if !in_block
            in_block = true
            block_start = i
        end
    end
    if in_block && last_above - block_start > min_size
        push!(blocks, shrink_block(block_start:last_above, 0.1, 0.9))
    end

    return blocks
end

function find_zeros(ts, data, trigger=0)
    crossings = Float64[]
    for i in 1:length(data) - 1
        d1 = data[i]
        d2 = data[i + 1]
        if d1 <= 0 && d2 > 0
            t1 = ts[i]
            t2 = ts[i + 1]
            push!(crossings, (t1 * d2 - t2 * d1) / (d2 - d1))
        end
    end
    return crossings
end

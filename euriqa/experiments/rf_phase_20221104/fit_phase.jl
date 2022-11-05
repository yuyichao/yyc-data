#!/usr/bin/julia

using HDF5

load_data(fname) = h5open(fname) do io
    return (time=read(io, "time"), value=read(io, "value"),
            tscale=read(io, "tscale"), vscale=read(io, "vscale"))
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
                    push!(blocks, (block_start, last_above))
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
        push!(blocks, (block_start, last_above))
    end

    return blocks
end

function shrink_block(block, r_start, r_end)
    i = first(block)
    l = last(block)
    new_i = round(Int, i * (1 - r_start) + l * r_start)
    new_l = round(Int, i * (1 - r_end) + l * r_end)
    return new_i:new_l
end

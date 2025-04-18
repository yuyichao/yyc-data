#!/usr/bin/julia

function redistricute_blocks!(blocks)
    block, idx = findmax(blocks)
    blocks[idx] = 0
    nblocks = length(blocks)
    for i in 1:block
        idx = idx % nblocks + 1
        blocks[idx] += 1
    end
end

function find_loop(file)
    blocks = parse.(Int, split(read(file, String)))
    seen = Set{Vector{Int}}()
    c = 0
    while true
        if blocks in seen
            return c
        end
        push!(seen, copy(blocks))
        redistricute_blocks!(blocks)
        c += 1
    end
end

@show find_loop(ARGS[1])

#!/usr/bin/julia

struct MoveInst
    from::Int
    to::Int
    count::Int
end

function move_crates!(stacks, inst)
    from = stacks[inst.from]
    to = stacks[inst.to]
    for i in 1:inst.count
        push!(to, pop!(from))
    end
end

function top_crate(file)
    lines = eachline(file)
    stacks = Vector{String}[]
    for line in lines
        if isempty(line)
            break
        end
        for m in eachmatch(r"\[([A-Z])\]", line)
            stacknum = (m.offset + 3) ÷ 4
            name = m[1]
            nstacks = length(stacks)
            if nstacks < stacknum
                resize!(stacks, stacknum)
                for i in nstacks + 1:stacknum
                    stacks[i] = String[]
                end
            end
            pushfirst!(stacks[stacknum], name)
        end
    end

    for line in lines
        m = match(r"move (\d+) from (\d+) to (\d+)", line)
        move_crates!(stacks, MoveInst(parse(Int, m[2]), parse(Int, m[3]),
                                      parse(Int, m[1])))
    end
    return join([stack[end] for stack in stacks])
end

@show top_crate(ARGS[1])

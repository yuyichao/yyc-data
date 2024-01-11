#!/usr/bin/julia

function sum_value(nums::AbstractVector{T} where T<:Integer)
    nchildren = popfirst!(nums)
    nmeta = popfirst!(nums)
    if nchildren == 0
        s = 0
        for i in 1:nmeta
            s += popfirst!(nums)
        end
        return s
    end
    s = 0
    children = [sum_value(nums) for i in 1:nchildren]
    for i in 1:nmeta
        idx = popfirst!(nums)
        if idx < 1 || idx > nchildren
            continue
        end
        s += children[idx]
    end
    return s
end

function sum_value(file)
    nums = parse.(Int, split(read(file, String)))
    s = sum_value(nums)
    @assert isempty(nums)
    return s
end

@show sum_value(ARGS[1])

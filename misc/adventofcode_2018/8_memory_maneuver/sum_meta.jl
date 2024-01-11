#!/usr/bin/julia

function sum_meta(nums::AbstractVector{T} where T<:Integer)
    nchildren = popfirst!(nums)
    nmeta = popfirst!(nums)
    s = 0
    for i in 1:nchildren
        s += sum_meta(nums)
    end
    for i in 1:nmeta
        s += popfirst!(nums)
    end
    return s
end

function sum_meta(file)
    nums = parse.(Int, split(read(file, String)))
    s = sum_meta(nums)
    @assert isempty(nums)
    return s
end

@show sum_meta(ARGS[1])

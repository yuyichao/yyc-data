#!/usr/bin/julia

function is_triangle(line)
    nums = sort!(parse.(Int, split(line)))
    @assert length(nums) == 3
    return nums[1] + nums[2] > nums[3]
end

@show count(is_triangle, eachline(ARGS[1]))

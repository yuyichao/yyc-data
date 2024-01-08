#!/usr/bin/julia

function calc_area(line)
    nums = parse.(Int, split(line, 'x'))
    sort!(nums)
    return 3 * nums[1] * nums[2] + 2 * nums[1] * nums[3] + 2 * nums[2] * nums[3]
end

@show sum(calc_area, eachline(ARGS[1]))

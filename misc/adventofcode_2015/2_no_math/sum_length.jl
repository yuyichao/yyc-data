#!/usr/bin/julia

function calc_length(line)
    nums = parse.(Int, split(line, 'x'))
    sort!(nums)
    return 2 * nums[1] + 2 * nums[2] + nums[1] * nums[2] * nums[3]
end

@show sum(calc_length, eachline(ARGS[1]))

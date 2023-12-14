#!/usr/bin/julia

function predict(nums)
    if length(nums) == 1
        return nums[1]
    end
    if all(==(0), nums)
        return 0
    end
    return nums[1] - predict(diff(nums))
end

function sum_prev(file)
    s = 0
    for line in eachline(file)
        s += predict(parse.(Int, split(line)))
    end
    return s
end

@show sum_prev(ARGS[1])

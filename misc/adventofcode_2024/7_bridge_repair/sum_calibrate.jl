#!/usr/bin/julia

function check_exp(tgt, nums, ops)
    n = length(nums)
    v = nums[1]
    for idx in 1:n - 1
        if v > tgt
            return false
        end
        v2 = nums[idx + 1]
        if ((ops >> (idx - 1)) & 1) != 0
            v *= v2
        else
            v += v2
        end
    end
    return tgt == v
end

function check_line(tgt, nums)
    n = length(nums)
    for i in 0:(1 << (n - 1)) - 1
        if check_exp(tgt, nums, i)
            return true
        end
    end
    return false
end

function sum_calibrate(file)
    s = 0
    for line in eachline(file)
        tgt_str, num_str = split(line, ":")
        tgt = parse(Int, tgt_str)
        nums = parse.(Int, split(num_str))
        if check_line(tgt, nums)
            s += tgt
        end
    end
    return s
end

@show sum_calibrate(ARGS[1])

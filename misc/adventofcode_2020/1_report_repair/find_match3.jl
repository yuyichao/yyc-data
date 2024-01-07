#!/usr/bin/julia

function find_match3(file)
    nums = [parse(Int, line) for line in eachline(file)]
    nnums = length(nums)
    for i in 1:nnums
        n1 = nums[i]
        for j in 1:nnums
            if j == i
                continue
            end
            n2 = nums[j]
            for k in 1:nnums
                if k == i || k == j
                    continue
                end
                n3 = nums[k]
                if n1 + n2 + n3 == 2020
                    return n1 * n2 * n3
                end
            end
        end
    end
end

@show find_match3(ARGS[1])

#!/usr/bin/julia

function find_divide(nums)
    n = length(nums)
    for i in 1:n
        n1 = nums[i]
        for j in i + 1:n
            n2 = nums[j]
            if n1 % n2 == 0
                return n1 ÷ n2
            elseif n2 % n1 == 0
                return n2 ÷ n1
            end
        end
    end
end

function sum_divide(file)
    return sum(find_divide(parse.(Int, split(line))) for line in eachline(file))
end

@show sum_divide(ARGS[1])

#!/usr/bin/julia

function is_triangle(nums)
    nums = sort(nums)
    return nums[1] + nums[2] > nums[3]
end

function count_triangles(file)
    nums = [parse.(Int, split(line)) for line in eachline(file)]
    c = 0
    for i in 1:3:length(nums)
        c += is_triangle((nums[i][1], nums[i + 1][1], nums[i + 2][1]))
        c += is_triangle((nums[i][2], nums[i + 1][2], nums[i + 2][2]))
        c += is_triangle((nums[i][3], nums[i + 1][3], nums[i + 2][3]))
    end
    return c
end

@show count_triangles(ARGS[1])

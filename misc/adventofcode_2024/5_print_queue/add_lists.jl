#!/usr/bin/julia

function add_list(file)
    lines = eachline(file)
    orders = Set{NTuple{2,Int}}()
    for line in lines
        if isempty(line)
            break
        end
        x, y = parse.(Int, split(line, "|"))
        push!(orders, (x, y))
    end

    s = 0
    for line in lines
        nums = parse.(Int, split(line, ","))
        for i in 2:length(nums)
            for j in 1:i - 1
                if (nums[i], nums[j]) in orders
                    @goto next_line
                end
            end
        end
        s += nums[end ÷ 2 + 1]
        @label next_line
    end
    return s
end

@show add_list(ARGS[1])

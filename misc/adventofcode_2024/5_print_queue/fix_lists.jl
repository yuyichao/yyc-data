#!/usr/bin/julia

function add_list(file)
    lines = eachline(file)
    orders = [Set{Int}() for i in 1:99]
    for line in lines
        if isempty(line)
            break
        end
        x, y = parse.(Int, split(line, "|"))
        push!(orders[y], x)
    end

    s = 0
    for line in lines
        nums = parse.(Int, split(line, ","))
        for i in 2:length(nums)
            for j in 1:i - 1
                if nums[i] in orders[nums[j]]
                    @goto fix_line
                end
            end
        end
        continue
        @label fix_line
        nums2 = Set(nums)
        empty!(nums)
        function add_num(n)
            @assert n in nums2
            delete!(nums2, n)
            for b in orders[n]
                if b in nums2
                    add_num(b)
                end
            end
            push!(nums, n)
        end
        while !isempty(nums2)
            add_num(first(nums2))
        end
        s += nums[end ÷ 2 + 1]
    end
    return s
end

@show add_list(ARGS[1])

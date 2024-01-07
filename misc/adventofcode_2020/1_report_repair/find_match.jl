#!/usr/bin/julia

function find_match(file)
    nums = Set(parse(Int, line) for line in eachline(file))
    for num in nums
        if (2020 - num) in nums
            return num * (2020 - num)
        end
    end
end

@show find_match(ARGS[1])

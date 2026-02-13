#!/usr/bin/julia

function sum_results(file)
    lines = [split(line) for line in eachline(file)]
    ops = lines[end]
    numbers = [begin
                   @assert length(line) == length(ops)
                   parse.(Int, line)
               end for line in @view lines[1:end - 1]]
    s = 0
    for i in 1:length(ops)
        if ops[i] == "+"
            s += sum(nums[i] for nums in numbers)
        else
            @assert ops[i] == "*"
            s += prod(nums[i] for nums in numbers)
        end
    end
    return s
end

@show sum_results(ARGS[1])

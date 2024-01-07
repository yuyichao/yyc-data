#!/usr/bin/julia

function compute_power(file)
    counts = Dict{Int,Vector{Int}}()
    for line in eachline(file)
        for (i, c) in enumerate(line)
            get!(()->zeros(Int, 2), counts, i)[c == '0' ? 1 : 2] += 1
        end
    end
    gamma_rate = String([counts[i][1] > counts[i][2] ? '0' : '1' for i in 1:length(counts)])
    epsilon_rate = String([counts[i][1] < counts[i][2] ? '0' : '1' for i in 1:length(counts)])
    return parse(Int, gamma_rate, base=2) * parse(Int, epsilon_rate, base=2)
end

@show compute_power(ARGS[1])

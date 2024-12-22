#!/usr/bin/julia

mix(s, v) = s ⊻ v
prune(s) = s % 16777216

function next_secret(s)
    s = prune(mix(s, s * 64))
    s = prune(mix(s, s ÷ 32))
    s = prune(mix(s, s * 2048))
    return s
end

function nth_secret(s, n)
    for i in 1:n
        s = next_secret(s)
    end
    return s
end

function sum_secret(file)
    s = 0
    for line in eachline(file)
        s += nth_secret(parse(Int, line), 2000)
    end
    return s
end

@show sum_secret(ARGS[1])

#!/usr/bin/julia

function decode_rep(file)
    res = Dict{Int,Dict{Char,Int}}()
    maxlen = 0
    for line in eachline(file)
        maxlen = max(length(line), maxlen)
        for (i, c) in enumerate(line)
            r = get!(Dict{Char,Int}, res, i)
            r[c] = get(r, c, 0) + 1
        end
    end
    chars = [findmin(res[i])[2] for i in 1:maxlen]
    return String(chars)
end

@show decode_rep(ARGS[1])

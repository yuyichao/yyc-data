#!/usr/bin/julia

const cache = Dict{String,Int}()

function count_pattern(pattern, flags, maxlen)
    return get!(cache, pattern) do
        l = length(pattern)
        s = 0
        for i in 1:min(l, maxlen)
            if pattern[1:i] in flags
                if i == l
                    s += 1
                else
                    s += count_pattern(pattern[i + 1:l], flags, maxlen)
                end
            end
        end
        return s
    end
end

function count_all_pattern(file, file_tgt)
    flags = Set(strip.(split(read(file, String), ",")))
    maxlen = maximum(length(f) for f in flags)
    s = 0
    for line in eachline(file_tgt)
        s += count_pattern(line, flags, maxlen)
    end
    return s
end

@show count_all_pattern(ARGS[1], ARGS[2])

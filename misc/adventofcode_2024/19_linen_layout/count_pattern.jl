#!/usr/bin/julia

function check_pattern(pattern, flags, maxlen)
    l = length(pattern)
    if l == 0
        return true
    end
    for i in 1:min(l, maxlen)
        if pattern[1:i] in flags && check_pattern(pattern[i + 1:l], flags, maxlen)
            return true
        end
    end
    return false
end

function count_pattern(file, file_tgt)
    flags = Set(strip.(split(read(file, String), ",")))
    maxlen = maximum(length(f) for f in flags)
    s = 0
    for line in eachline(file_tgt)
        if check_pattern(line, flags, maxlen)
            s += 1
        end
    end
    return s
end

@show count_pattern(ARGS[1], ARGS[2])

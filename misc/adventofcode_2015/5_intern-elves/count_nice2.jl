#!/usr/bin/julia

function is_nice(s)
    return match(r"(..).*\1", s) !== nothing && match(r"(.).\1", s) !== nothing
end

function count_nice(file)
    return count(is_nice, eachline(file))
end

@show count_nice(ARGS[1])

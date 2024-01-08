#!/usr/bin/julia

function sum_digits(file)
    str = strip(read(file, String))
    prev = str[end]
    s = 0
    for c in str
        if prev == c
            s += c - '0'
        end
        prev = c
    end
    return s
end

@show sum_digits(ARGS[1])

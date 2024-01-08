#!/usr/bin/julia

function sum_digits(file)
    str = strip(read(file, String))
    len = length(str)
    skip = len ÷ 2
    s = 0
    for (i, c) in enumerate(str)
        i2 = (i + skip - 1) % len + 1
        if c == str[i2]
            s += c - '0'
        end
    end
    return s
end

@show sum_digits(ARGS[1])

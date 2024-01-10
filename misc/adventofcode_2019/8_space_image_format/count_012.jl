#!/usr/bin/julia

function count_012(file)
    str = strip(read(file, String))
    len = length(str)
    min_c0 = 25 * 6
    min_c1 = 0
    min_c2 = 0
    for i in 1:(25 * 6):len
        substr = str[i:i + 25 * 6 - 1]
        c0 = count(==('0'), substr)
        if c0 < min_c0
            min_c0 = c0
            min_c1 = count(==('1'), substr)
            min_c2 = count(==('2'), substr)
        end
    end
    return min_c1 * min_c2
end

@show count_012(ARGS[1])

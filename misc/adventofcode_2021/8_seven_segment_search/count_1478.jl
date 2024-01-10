#!/usr/bin/julia

function count_1478(file)
    c = 0
    for line in eachline(file)
        _, line = split(line, '|')
        for digit in split(line)
            n = length(digit)
            c += (n == 2 || n == 4 || n == 3 || n == 7)
        end
    end
    return c
end

@show count_1478(ARGS[1])

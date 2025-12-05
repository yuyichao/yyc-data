#!/usr/bin/julia

function count_0(file)
    pos = 50
    count = 0
    for line in eachline(file)
        sign = line[1] == 'R' ? 1 : -1
        step = parse(Int, line[2:end])
        for i in 1:step
            pos = (pos + sign) % 100
            if pos == 0
                count += 1
            end
        end
    end
    return count
end

@show count_0(ARGS[1])

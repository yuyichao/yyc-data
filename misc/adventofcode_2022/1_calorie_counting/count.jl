#!/usr/bin/julia

function count(file)
    m = 0
    c = 0
    for line in eachline(file)
        if isempty(line)
            m = max(m, c)
            c = 0
        else
            c += parse(Int, line)
        end
    end
    return max(m, c)
end

@show count(ARGS[1])

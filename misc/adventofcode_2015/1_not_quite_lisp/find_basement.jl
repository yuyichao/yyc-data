#!/usr/bin/julia

function count_floor(file)
    f = 0
    for (i, c) in enumerate(strip(read(file, String)))
        if c == '('
            f += 1
        else
            f -= 1
        end
        if f == -1
            return i
        end
    end
end

@show count_floor(ARGS[1])

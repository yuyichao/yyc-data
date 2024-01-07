#!/usr/bin/julia

function compute_position(file)
    h = 0
    d = 0
    aim = 0
    for line in eachline(file)
        m = match(r"(forward|down|up) (\d+)", line)
        s = parse(Int, m[2])
        if m[1] == "forward"
            h += s
            d += aim * s
        elseif m[1] == "down"
            aim += s
        elseif m[1] == "up"
            aim -= s
        else
            error()
        end
    end
    return h * d
end

@show compute_position(ARGS[1])

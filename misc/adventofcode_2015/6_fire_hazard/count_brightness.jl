#!/usr/bin/julia

function count_brightness(file)
    lights = zeros(Int, 1000, 1000)
    for line in eachline(file)
        m = match(r"(toggle|turn off|turn on) (\d+),(\d+) through (\d+),(\d+)", line)
        x1 = parse(Int, m[2]) + 1
        y1 = parse(Int, m[3]) + 1
        x2 = parse(Int, m[4]) + 1
        y2 = parse(Int, m[5]) + 1

        x1, x2 = min(x1, x2), max(x1, x2)
        y1, y2 = min(y1, y2), max(y1, y2)

        if m[1] == "toggle"
            for x in x1:x2
                for y in y1:y2
                    lights[y, x] += 2
                end
            end
        elseif m[1] == "turn off"
            for x in x1:x2
                for y in y1:y2
                    lights[y, x] = max(lights[y, x] - 1, 0)
                end
            end
        elseif m[1] == "turn on"
            for x in x1:x2
                for y in y1:y2
                    lights[y, x] += 1
                end
            end
        end
    end
    return sum(lights)
end

@show count_brightness(ARGS[1])

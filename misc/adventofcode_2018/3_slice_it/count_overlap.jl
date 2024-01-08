#!/usr/bin/julia

function count_overlap(file)
    M = zeros(Int, 1000, 1000)
    for line in eachline(file)
        x0, y0, dx, dy = parse.(Int, match(r"#\d+ @ (\d+),(\d+): (\d+)x(\d+)", line))
        for x in x0 + 1:x0 + dx
            for y in y0 + 1:y0 + dy
                M[y, x] += 1
            end
        end
    end
    return count(>(1), M)
end

@show count_overlap(ARGS[1])

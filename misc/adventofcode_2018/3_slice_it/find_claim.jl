#!/usr/bin/julia

function find_claim(file)
    M = zeros(Int, 1000, 1000)
    for line in eachline(file)
        x0, y0, dx, dy = parse.(Int, match(r"#\d+ @ (\d+),(\d+): (\d+)x(\d+)", line))
        for x in x0 + 1:x0 + dx
            for y in y0 + 1:y0 + dy
                M[y, x] += 1
            end
        end
    end
    for line in eachline(file)
        id, x0, y0, dx, dy = parse.(Int, match(r"#(\d+) @ (\d+),(\d+): (\d+)x(\d+)", line))
        has_overlap = false
        for x in x0 + 1:x0 + dx
            for y in y0 + 1:y0 + dy
                if M[y, x] > 1
                    has_overlap = true
                end
            end
        end
        if !has_overlap
            return id
        end
    end
end

@show find_claim(ARGS[1])

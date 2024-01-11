#!/usr/bin/julia

function count_pixels(file)
    screen = zeros(Bool, 6, 50)
    for line in eachline(file)
        m = match(r"rect (\d+)x(\d+)", line)
        if m !== nothing
            screen[1:parse(Int, m[2]), 1:parse(Int, m[1])] .= true
            continue
        end
        m = match(r"rotate (row y|column x)=(\d+) by (\d+)", line)
        v1 = parse(Int, m[2]) + 1
        shift = parse(Int, m[3])
        if m[1] == "row y"
            shift = shift % size(screen, 2)
            if shift == 0
                continue
            end
            screen[v1, :] = [screen[v1, end - shift + 1:end]; screen[v1, 1:end - shift]]
        else
            @assert m[1] == "column x"
            shift = shift % size(screen, 1)
            if shift == 0
                continue
            end
            screen[:, v1] = [screen[end - shift + 1:end, v1]; screen[1:end - shift, v1]]
        end
    end
    return count(screen)
end

@show count_pixels(ARGS[1])

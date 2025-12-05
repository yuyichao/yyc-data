#!/usr/bin/julia

function count_paper(file)
    lines = readlines(file)
    sz = length(lines), length(lines[1])
    m = Int[try
                lines[i][j] == '@'
            catch
                0
            end
            for i in 0:length(lines) + 1, j in 0:length(lines[1]) + 1]
    c = 0
    for i in 2:size(m, 1) - 1
        for j in 2:size(m, 2) - 1
            if m[i, j] == 0
                continue
            end
            if m[i - 1, j - 1] + m[i - 1, j] + m[i - 1, j + 1] + m[i, j - 1] + m[i, j + 1] + m[i + 1, j - 1] + m[i + 1, j] + m[i + 1, j + 1] < 4
                c += 1
            end
        end
    end
    return c
end

@show count_paper(ARGS[1])

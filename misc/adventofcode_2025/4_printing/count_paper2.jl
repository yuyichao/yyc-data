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

    orig = sum(m)
    while true
        changed = false
        for i in 2:size(m, 1) - 1
            for j in 2:size(m, 2) - 1
                if m[i, j] == 0
                    continue
                end
                if m[i - 1, j - 1] + m[i - 1, j] + m[i - 1, j + 1] + m[i, j - 1] + m[i, j + 1] + m[i + 1, j - 1] + m[i + 1, j] + m[i + 1, j + 1] < 4
                    m[i, j] = 0
                    changed = true
                end
            end
        end
        if !changed
            break
        end
    end
    return orig - sum(m)
end

@show count_paper(ARGS[1])

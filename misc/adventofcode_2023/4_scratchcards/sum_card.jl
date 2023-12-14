#!/usr/bin/julia

function count_matching_line(line)
    _, results = split(line, ":")
    win, mine = split(results, "|")
    win = Set(split(win))
    n = 0
    for v in split(mine)
        if v in win
            n += 1
        end
    end
    return n
end

function count_card(file)
    matching = [count_matching_line(line) for line in eachline(file)]
    nlines = length(matching)
    copies = ones(Int, nlines)
    for i in 1:nlines
        nm = matching[i]
        i1 = i + 1
        i2 = min(i + nm, nlines)
        if i2 < i1
            continue
        end
        copies[i1:i2] .+= copies[i]
    end
    return sum(copies)
end

@show count_card(ARGS[1])

#!/usr/bin/julia

function count_answers(file)
    c = 0
    set = Set{Char}()
    for line in eachline(file)
        if isempty(line)
            c += length(set)
            empty!(set)
        else
            union!(set, line)
        end
    end
    c += length(set)
    return c
end

@show count_answers(ARGS[1])

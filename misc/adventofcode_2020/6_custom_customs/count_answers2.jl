#!/usr/bin/julia

function count_answers(file)
    c = 0
    is_first = true
    set = Set{Char}()
    for line in eachline(file)
        if isempty(line)
            c += length(set)
            empty!(set)
            is_first = true
        elseif is_first
            union!(set, line)
            is_first = false
        else
            intersect!(set, line)
        end
    end
    c += length(set)
    return c
end

@show count_answers(ARGS[1])

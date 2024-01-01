#!/usr/bin/julia

function is_nice(s)
    last_c = '\0'
    vowels_count = 0
    has_double = false
    for c in s
        if c in "aeiou"
            vowels_count += 1
        end
        if last_c == c
            has_double = true
        end
        if (last_c, c) in (('a', 'b'), ('c', 'd'), ('p', 'q'), ('x', 'y'))
            return false
        end
        last_c = c
    end
    return has_double && vowels_count >= 3
end

function count_nice(file)
    return count(is_nice, eachline(file))
end

@show count_nice(ARGS[1])

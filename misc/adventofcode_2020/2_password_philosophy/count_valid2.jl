#!/usr/bin/julia

function check_password(pw, lo, hi, c0)
    return (pw[lo] == c0) ⊻ (pw[hi] == c0)
end

function count_valid(file)
    c = 0
    for line in eachline(file)
        m = match(r"(\d+)-(\d+) (.): *(.+)", line)
        if check_password(m[4], parse(Int, m[1]), parse(Int, m[2]), m[3][1])
            c += 1
        end
    end
    return c
end

@show count_valid(ARGS[1])

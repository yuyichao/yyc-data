#!/usr/bin/julia

function is_overlap(rng1, rng2)
    return !(rng1[2] < rng2[1] || rng2[2] < rng1[1])
end

function count_overlap(file)
    return count(eachline(file)) do line
        m = match(r"(\d+)-(\d+),(\d+)-(\d+)", line)
        return is_overlap((parse(Int, m[1]), parse(Int, m[2])),
                          (parse(Int, m[3]), parse(Int, m[4])))
    end
end

@show count_overlap(ARGS[1])

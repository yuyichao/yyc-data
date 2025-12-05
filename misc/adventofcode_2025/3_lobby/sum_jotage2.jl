#!/usr/bin/julia

function find_max_digit(ds, rng)
    maxd = 0
    maxi = 0
    for i in rng
        d = ds[i]
        if d > maxd
            maxi = i
            maxd = d
            if d == 9
                break
            end
        end
    end
    return maxi, maxd
end

function find_max_jotage(ds, nc)
    nd = length(ds)
    v = 0
    i0 = 1
    for p in 1:nc
        maxi, maxd = find_max_digit(ds, i0:nd - nc + p)
        i0 = maxi + 1
        v = v * 10 + maxd
    end
    return v
end

function sum_jotage(file)
    s = 0
    for line in eachline(file)
        s += find_max_jotage([c - '0' for c in line], 12)
    end
    return s
end

@show sum_jotage(ARGS[1])

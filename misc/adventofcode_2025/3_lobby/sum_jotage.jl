#!/usr/bin/julia

function sum_jotage(file)
    s = 0
    for line in eachline(file)
        ds = [c - '0' for c in line]
        nd = length(ds)
        maxd1 = 0
        maxi1 = 0
        for i in 1:nd - 1
            d = ds[i]
            if d > maxd1
                maxi1 = i
                maxd1 = d
                if d == 9
                    break
                end
            end
        end
        maxd2 = 0
        for i in maxi1 + 1:nd
            d = ds[i]
            if d > maxd2
                maxd2 = d
                if d == 9
                    break
                end
            end
        end
        s += maxd1 * 10 + maxd2
    end
    return s
end

@show sum_jotage(ARGS[1])

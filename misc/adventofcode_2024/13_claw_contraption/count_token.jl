#!/usr/bin/julia

function check_game(a, b, p)
    pa = (0, 0)
    ta = 0
    mint = -1
    while pa[1] <= p[1] && pa[2] <= p[2]
        pb = pa
        tb = ta
        while pb[1] <= p[1] && pb[2] <= p[2]
            if pb == p
                if mint == -1 || mint > tb
                    mint = tb
                    break
                end
            end
            pb = pb .+ b
            tb += 1
        end
        pa = pa .+ a
        ta += 3
    end
    return mint == -1 ? 0 : mint
end

function count_token(file)
    lines = readlines(file)
    s = 0
    for i in 1:4:length(lines)
        ma = match(r"Button A: X\+(\d*), Y\+(\d*)", lines[i])
        mb = match(r"Button B: X\+(\d*), Y\+(\d*)", lines[i + 1])
        mp = match(r"Prize: X=(\d*), Y=(\d*)", lines[i + 2])
        a = (parse(Int, ma[1]), parse(Int, ma[2]))
        b = (parse(Int, mb[1]), parse(Int, mb[2]))
        p = (parse(Int, mp[1]), parse(Int, mp[2]))
        s += check_game(a, b, p)
    end
    return s
end

@show count_token(ARGS[1])

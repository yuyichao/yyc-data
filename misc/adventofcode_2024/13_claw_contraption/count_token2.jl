#!/usr/bin/julia

function check_game(a, b, p)
    s = [a[1] b[1]
         a[2] b[2]] \ [p[1], p[2]]
    s = (round(Int, s[1]), round(Int, s[2]))
    if (a .* s[1] .+ b .* s[2]) == p
        return s[1] * 3 + s[2]
    end
    return 0
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
        @assert !(a[1] / b[1] ≈ a[2] / b[2])
        # p = (parse(Int, mp[1]) + 10000000000000, parse(Int, mp[2]) + 10000000000000)
        p = (parse(Int, mp[1]) + 10000000000000, parse(Int, mp[2]) + 10000000000000)
        s += check_game(a, b, p)
    end
    return s
end

@show count_token(ARGS[1])

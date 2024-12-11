#!/usr/bin/julia

function evolve!(tgt, src)
    empty!(tgt)
    for n in src
        if n == 0
            push!(tgt, 1)
            continue
        end
        s = string(n)
        l = length(s)
        if l % 2 == 0
            push!(tgt, parse(BigInt, s[1:l ÷ 2]))
            push!(tgt, parse(BigInt, s[l ÷ 2 + 1:l]))
        else
            push!(tgt, n * 2024)
        end
    end
end

function evolve_n!(tgt, src, n)
    for i in 1:n
        evolve!(tgt, src)
        tgt, src = src, tgt
    end
    return src
end

@show length(evolve_n!(BigInt[], BigInt[5, 62914, 65, 972, 0, 805922, 6521, 1639064], 25))

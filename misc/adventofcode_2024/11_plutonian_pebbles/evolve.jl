#!/usr/bin/julia

function evolve!(tgt, src)
    empty!(tgt)
    for n in src
        if n == 0
            push!(tgt, 1)
            continue
        end
        if 10 <= n < 100
            push!(tgt, n ÷ 10)
            push!(tgt, n % 10)
        elseif 1000 <= n < 10000
            push!(tgt, n ÷ 100)
            push!(tgt, n % 100)
        elseif 100000 <= n < 1000000
            push!(tgt, n ÷ 1000)
            push!(tgt, n % 1000)
        elseif 10000000 <= n < 100000000
            push!(tgt, n ÷ 10000)
            push!(tgt, n % 10000)
        elseif 1000000000 <= n < 10000000000
            push!(tgt, n ÷ 100000)
            push!(tgt, n % 100000)
        elseif 100000000000 <= n < 1000000000000
            push!(tgt, n ÷ 1000000)
            push!(tgt, n % 1000000)
        else
            @assert n < 1000000000000
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

@show length(evolve_n!(Int[], Int[5, 62914, 65, 972, 0, 805922, 6521, 1639064], 25))

#!/usr/bin/julia

const evolve_rule = Dict{Int,Vector{Int}}()

function collect_rule(n)
    get!(evolve_rule, n) do
        if n == 0
            return [(1)]
        end
        if 10 <= n < 100
            return [(n ÷ 10), (n % 10)]
        elseif 1000 <= n < 10000
            return [(n ÷ 100), (n % 100)]
        elseif 100000 <= n < 1000000
            return [(n ÷ 1000), (n % 1000)]
        elseif 10000000 <= n < 100000000
            return [(n ÷ 10000), (n % 10000)]
        elseif 1000000000 <= n < 10000000000
            return [(n ÷ 100000), (n % 100000)]
        elseif 100000000000 <= n < 1000000000000
            return [(n ÷ 1000000), (n % 1000000)]
        else
            @assert n < 1000000000000
            return [(n * 2024)]
        end
    end
    return n
end

const input = [5, 62914, 65, 972, 0, 805922, 6521, 1639064]
for n in input
    collect_rule(n)
end

while true
    nrules = length(evolve_rule)
    for (k, v) in evolve_rule
        for n in v
            collect_rule(n)
        end
    end
    if length(evolve_rule) == nrules
        break
    end
end

function evolve!(tgt, src)
    empty!(tgt)
    for (n, r) in src
        for o in evolve_rule[n]
            tgt[o] = get(tgt, o, 0) + r
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

@show sum(r for (n, r) in evolve_n!(Dict{Int,Int}(),
                                    Dict([n=>1 for n in input]), 75))

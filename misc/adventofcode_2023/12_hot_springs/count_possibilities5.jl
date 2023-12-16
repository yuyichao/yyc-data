#!/usr/bin/julia

using Combinatorics

# 1: operational, 0: not
@inline function count_startat(condition, groups, ci, gi, ncond, ngroup,
                               rest_count, cache)
    g = @inbounds groups[gi]
    @inbounds for i in 1:g
        idx = ci + i - 1
        if idx > ncond || condition[idx] == 1
            return 0
        end
    end
    idx = ci + g
    @inbounds if idx <= ncond && condition[idx] == 0
        return 0
    end
    return count_findstart(condition, groups, idx + 1, gi + 1, ncond,
                           ngroup, rest_count, cache)
end

@inline function count_findstart(condition, groups, ci, gi, ncond, ngroup,
                                 rest_count, cache)
    return get!(cache, (ci, gi)) do
        if gi > ngroup
            @inbounds for idx in ci:ncond
                if condition[idx] == 0
                    return 0
                end
            end
            return 1
        end
        count = 0
        @inbounds for idx in ci:(ncond - rest_count[gi] + 1)
            c = condition[idx]
            if c == 1
                continue
            elseif c == 0
                return count + count_startat(condition, groups, idx, gi,
                                             ncond, ngroup, rest_count, cache)
            else
                count += count_startat(condition, groups, idx, gi,
                                       ncond, ngroup, rest_count, cache)
            end
        end
        return count
    end
end

function enumerate_conditions(condition, groups)
    condition = repeat([condition; -1], 5)[1:end - 1]
    groups = repeat(groups, 5)
    rest_count = copy(groups)
    ngroup = length(groups)
    for i in ngroup - 1:-1:1
        rest_count[i] += rest_count[i + 1] + 1
    end
    return count_findstart(condition, groups, 1, 1, length(condition),
                           ngroup, rest_count, Dict{NTuple{2,Int},Int}())
end

function count_possibilities(file)
    s = 0
    i = 0
    for line in eachline(file)
        cond, group = split(line)
        i += 1
        s += enumerate_conditions([c == '.' ? 1 : (c == '#' ? 0 : -1) for c in cond],
                                  parse.(Int, split(group, ',')))
    end
    return s
end

@show count_possibilities(ARGS[1])

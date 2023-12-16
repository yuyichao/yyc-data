#!/usr/bin/julia

using Combinatorics

# 1 operational, 0 not
function count_groups(condition)
    groups = Int[]
    state = 1
    first_broken = 0
    for i in 1:length(condition)
        c = condition[i]
        if c == 1
            if state == 0
                state = 1
                push!(groups, i - first_broken)
            end
        elseif state == 1
            state = 0
            first_broken = i
        end
    end
    if state == 0
        push!(groups, length(condition) + 1 - first_broken)
    end
    return groups
end

function enumerate_conditions(condition, groups)
    remain_broken = sum(groups)
    blank = copy(condition)
    buff = copy(condition)
    unknown_pos = Int[]
    for i in 1:length(condition)
        c = condition[i]
        if c == 0
            remain_broken -= 1
        elseif c == -1
            blank[i] = 1
            push!(unknown_pos, i)
        end
    end

    s = 0
    for extra_broken in combinations(unknown_pos, remain_broken)
        buff .= blank
        for i in extra_broken
            buff[i] = 0
        end
        if count_groups(buff) == groups
            s += 1
        end
    end
    return s
end

function count_possibilities(file)
    s = 0
    for line in eachline(file)
        cond, group = split(line)
        s += enumerate_conditions([c == '.' ? 1 : (c == '#' ? 0 : -1) for c in cond],
                                  parse.(Int, split(group, ',')))
    end
    return s
end

@show count_possibilities(ARGS[1])

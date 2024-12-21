#!/usr/bin/julia

const num_pos = Dict('7'=>(1, 1), '8'=>(1, 2), '9'=>(1, 3),
                     '4'=>(2, 1), '5'=>(2, 2), '6'=>(2, 3),
                     '1'=>(3, 1), '2'=>(3, 2), '3'=>(3, 3),
                     '0'=>(4, 2), 'A'=>(4, 3))
const dir_pos = Dict('^'=>(1, 2), 'A'=>(1, 3),
                     '<'=>(2, 1), 'v'=>(2, 2), '>'=>(2, 3))

function get_keypress_options(str, pos_map, invalid_pos)
    cur = pos_map['A']
    options = Vector{String}[]
    for c in str
        next = pos_map[c]
        if cur[1] == next[1]
            if cur[2] == next[2]
                push!(options, ["A"])
            elseif cur[2] < next[2]
                push!(options, [">"^(next[2] - cur[2]) * "A"])
            else
                @assert cur[2] > next[2]
                push!(options, ["<"^(cur[2] - next[2]) * "A"])
            end
        elseif cur[2] == next[2]
            if cur[1] < next[1]
                push!(options, ["v"^(next[1] - cur[1]) * "A"])
            else
                @assert cur[1] > next[1]
                push!(options, ["^"^(cur[1] - next[1]) * "A"])
            end
        else
            if cur[1] < next[1]
                vs = "v"^(next[1] - cur[1])
            else
                @assert cur[1] > next[1]
                vs = "^"^(cur[1] - next[1])
            end
            if cur[2] < next[2]
                hs = ">"^(next[2] - cur[2])
            else
                @assert cur[2] > next[2]
                hs = "<"^(cur[2] - next[2])
            end
            if (cur[1], next[2]) == invalid_pos
                push!(options, [vs * hs * "A"])
            elseif (next[1], cur[2]) == invalid_pos
                push!(options, [hs * vs * "A"])
            else
                push!(options, [vs * hs * "A", hs * vs * "A"])
            end
        end
        cur = next
    end
    res = String[]
    min_len = typemax(Int)
    for strs in Iterators.product(options...)
        s = join(strs, "")
        min_len = min(min_len, length(s))
        push!(res, s)
    end
    return [s for s in res if length(s) == min_len]
end

function get_options(str)
    s = typemax(Int)
    for num1 in get_keypress_options(str, num_pos, (4, 1))
        for num2 in get_keypress_options(num1, dir_pos, (1, 1))
            for num3 in get_keypress_options(num2, dir_pos, (1, 1))
                s = min(s, length(num3))
            end
        end
    end
    return s
end

function sum_keys(file)
    s = 0
    for str in eachline(file)
        s += get_options(str) * parse(Int, str[1:end - 1])
    end
    return s
end

@show sum_keys(ARGS[1])

#!/usr/bin/julia

const num_pos = Dict('7'=>(1, 1), '8'=>(1, 2), '9'=>(1, 3),
                     '4'=>(2, 1), '5'=>(2, 2), '6'=>(2, 3),
                     '1'=>(3, 1), '2'=>(3, 2), '3'=>(3, 3),
                     '0'=>(4, 2), 'A'=>(4, 3))
const dir_pos = Dict('^'=>(1, 2), 'A'=>(1, 3),
                     '<'=>(2, 1), 'v'=>(2, 2), '>'=>(2, 3))

function get_keypress_options!(res, str, pos_map, invalid_pos, weight=1)
    cur = pos_map['A']
    function add_str(s)
        res[s] = get(res, s, 0) + weight
    end
    for c in str
        next = pos_map[c]
        if cur[1] == next[1]
            if cur[2] == next[2]
                add_str("A")
            elseif cur[2] < next[2]
                add_str(">"^(next[2] - cur[2]) * "A")
            else
                @assert cur[2] > next[2]
                add_str("<"^(cur[2] - next[2]) * "A")
            end
        elseif cur[2] == next[2]
            if cur[1] < next[1]
                add_str("v"^(next[1] - cur[1]) * "A")
            else
                @assert cur[1] > next[1]
                add_str("^"^(cur[1] - next[1]) * "A")
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
                add_str(vs * hs * "A")
            elseif (next[1], cur[2]) == invalid_pos
                add_str(hs * vs * "A")
            elseif cur[1] < next[1] && cur[2] < next[2]
                add_str(vs * hs * "A")
            elseif cur[1] < next[1] && cur[2] > next[2]
                add_str(hs * vs * "A")
            elseif cur[1] > next[1] && cur[2] < next[2]
                add_str(vs * hs * "A")
            else
                add_str(hs * vs * "A")
            end
        end
        cur = next
    end
    return res
end

function get_keypress_options2(prev)
    res = Dict{String,Int}()
    for (k, v) in prev
        get_keypress_options!(res, k, dir_pos, (1, 1), v)
    end
    return res
end

function get_options(str)
    num = get_keypress_options!(Dict{String,Int}(), str, num_pos, (4, 1))
    for i in 1:25
        num = get_keypress_options2(num)
    end
    return sum(length(k) * v for (k, v) in num)
end

function sum_keys(file)
    s = 0
    for str in eachline(file)
        s += get_options(str) * parse(Int, str[1:end - 1])
    end
    return s
end

@show sum_keys(ARGS[1])

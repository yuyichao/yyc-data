#!/usr/bin/julia

function run_intcode!(code)
    idx = 1
    while true
        c = code[idx]
        if c == 1
            p1 = code[idx + 1] + 1
            p2 = code[idx + 2] + 1
            p3 = code[idx + 3] + 1
            code[p3] = code[p1] + code[p2]
            idx += 4
        elseif c == 2
            p1 = code[idx + 1] + 1
            p2 = code[idx + 2] + 1
            p3 = code[idx + 3] + 1
            code[p3] = code[p1] * code[p2]
            idx += 4
        else
            @assert c == 99
            return
        end
    end
end

function find_result(orig_code, res)
    for noun in 0:99
        for verb in 0:99
            code = copy(orig_code)
            code[2] = noun
            code[3] = verb
            run_intcode!(code)
            if code[1] == res
                return noun * 100 + verb
            end
        end
    end
end

function run_intcode(file)
    code = parse.(Int, split(read(file, String), ','))
    return find_result(code, 19690720)
end

@show run_intcode(ARGS[1])

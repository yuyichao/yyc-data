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

function run_intcode(file)
    code = parse.(Int, split(read(file, String), ','))
    code[2] = 12
    code[3] = 2
    run_intcode!(code)
    return code[1]
end

@show run_intcode(ARGS[1])

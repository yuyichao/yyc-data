#!/usr/bin/julia

using Combinatorics

get_mode(c, i) = (c ÷ 10^(i + 1)) % 10

function get_param(code, idx, mode)
    c = code[idx]
    if mode == 0
        return code[c + 1]
    else
        @assert mode == 1
        return c
    end
end

function set_param(code, idx, mode, val)
    @assert mode == 0
    c = code[idx]
    code[c + 1] = val
    return
end

function run_intcode!(code, inputs, outputs)
    idx = 1
    while true
        c = code[idx]
        op = c % 100
        if op == 1
            set_param(code, idx + 3, get_mode(c, 3),
                      get_param(code, idx + 1, get_mode(c, 1)) +
                          get_param(code, idx + 2, get_mode(c, 2)))
            idx += 4
        elseif op == 2
            set_param(code, idx + 3, get_mode(c, 3),
                      get_param(code, idx + 1, get_mode(c, 1)) *
                          get_param(code, idx + 2, get_mode(c, 2)))
            idx += 4
        elseif op == 3
            set_param(code, idx + 1, get_mode(c, 1), take!(inputs))
            idx += 2
        elseif op == 4
            push!(outputs, get_param(code, idx + 1, get_mode(c, 1)))
            idx += 2
        elseif op == 5
            if get_param(code, idx + 1, get_mode(c, 1)) != 0
                idx = get_param(code, idx + 2, get_mode(c, 2)) + 1
            else
                idx += 3
            end
        elseif op == 6
            if get_param(code, idx + 1, get_mode(c, 1)) == 0
                idx = get_param(code, idx + 2, get_mode(c, 2)) + 1
            else
                idx += 3
            end
        elseif op == 7
            set_param(code, idx + 3, get_mode(c, 3),
                      get_param(code, idx + 1, get_mode(c, 1)) <
                          get_param(code, idx + 2, get_mode(c, 2)))
            idx += 4
        elseif op == 8
            set_param(code, idx + 3, get_mode(c, 3),
                      get_param(code, idx + 1, get_mode(c, 1)) ==
                          get_param(code, idx + 2, get_mode(c, 2)))
            idx += 4
        else
            @assert c == 99
            return
        end
    end
end

function run_amplifiers(code, phases)
    inputs = [Channel{Int}(3) for i in 1:5]
    tasks = Task[]

    for i in 1:5
        input = inputs[i]
        output = inputs[i % 5 + 1]
        push!(input, phases[i])
        if i == 1
            push!(input, 0)
        end
    end
    for i in 1:5
        let input = inputs[i], output = inputs[i % 5 + 1]
            t = Task(()->run_intcode!(copy(code), input, output))
            schedule(t)
            push!(tasks, t)
        end
    end
    for i in 1:5
        wait(tasks[i])
    end
    return take!(inputs[1])
end

function max_signal(file)
    code = parse.(Int, split(read(file, String), ','))
    return maximum(run_amplifiers(code, phases)
                   for phases in permutations(5:9))
end

@show max_signal(ARGS[1])

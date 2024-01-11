#!/usr/bin/julia

function max_registers(file)
    regs = Dict{String,Int}()
    for line in eachline(file)
        m = match(r"([a-z]+) (inc|dec) ([-+]?\d+) if ([a-z]+) (>|>=|<|<=|==|!=) ([-+]?\d+)", line)
        cond_v1 = get(regs, m[4], 0)
        cond_v2 = parse(Int, m[6])
        if m[5] == ">"
            cond = cond_v1 > cond_v2
        elseif m[5] == ">="
            cond = cond_v1 >= cond_v2
        elseif m[5] == "<"
            cond = cond_v1 < cond_v2
        elseif m[5] == "<="
            cond = cond_v1 <= cond_v2
        elseif m[5] == "=="
            cond = cond_v1 == cond_v2
        elseif m[5] == "!="
            cond = cond_v1 != cond_v2
        else
            error()
        end
        if !cond
            continue
        end

        v = parse(Int, m[3])
        if m[2] == "inc"
            regs[m[1]] = get(regs, m[1], 0) + v
        elseif m[2] == "dec"
            regs[m[1]] = get(regs, m[1], 0) - v
        else
            error()
        end
    end
    return findmax(regs)
end

@show max_registers(ARGS[1])

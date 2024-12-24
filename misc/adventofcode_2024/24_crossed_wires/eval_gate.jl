#!/usr/bin/julia

using Printf

struct Gate
    ty::Symbol
    in1::String
    in2::String
end

function eval_value(values, gates, name)
    if haskey(values, name)
        return values[name]
    end
    gate = gates[name]
    in1 = eval_value(values, gates, gate.in1)
    in2 = eval_value(values, gates, gate.in2)
    if gate.ty == :AND
        res = in1 & in2
    elseif gate.ty == :OR
        res = in1 | in2
    elseif gate.ty == :XOR
        res = in1 ⊻ in2
    else
        error("...")
    end
    values[name] = res
    return res
end

function eval_gate(file, file_init)
    values = Dict{String,Bool}()
    gates = Dict{String,Gate}()
    for line in eachline(file_init)
        name, val = split(line, ":")
        values[name] = parse(Int, val) != 0
    end
    for line in eachline(file)
        m = match(r"([^ ]*) (AND|OR|XOR) ([^ ]*) -> ([^ ]*)", line)
        gates[m[4]] = Gate(Symbol(m[2]), m[1], m[3])
    end
    zvs = Bool[]
    zv = 0
    while true
        zn = @sprintf("z%02d", zv)
        zv += 1
        if !haskey(gates, zn)
            break
        end
        push!(zvs, eval_value(values, gates, zn))
    end
    return reverse(zvs)
end

@show eval_gate(ARGS[1], ARGS[2])

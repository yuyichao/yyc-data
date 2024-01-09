#!/usr/bin/julia

mutable struct WireInfo
    w
    value_set::Bool
    value::UInt16
    WireInfo(w) = new(w, false, 0x0)
end

function get_value(wires, name::AbstractString)
    wire = wires[name]
    if !wire.value_set
        wire.value_set = true
        wire.value = eval(wire.w, wires)::UInt16
    end
    return wire.value
end

function get_value(wires, value::Integer)
    return UInt16(value)
end

struct Value
    x::Union{String,UInt16}
end
eval(w::Value, wires) = get_value(wires, w.x)

struct NotGate
    w::Union{String,UInt16}
end
eval(gate::NotGate, wires) = ~get_value(wires, gate.w)

struct AndGate
    w1::Union{String,UInt16}
    w2::Union{String,UInt16}
end
eval(gate::AndGate, wires) = get_value(wires, gate.w1) & get_value(wires, gate.w2)

struct OrGate
    w1::Union{String,UInt16}
    w2::Union{String,UInt16}
end
eval(gate::OrGate, wires) = get_value(wires, gate.w1) | get_value(wires, gate.w2)

struct LShiftGate
    w::Union{String,UInt16}
    v::Union{String,UInt16}
end
eval(gate::LShiftGate, wires) = get_value(wires, gate.w) << get_value(wires, gate.v)

struct RShiftGate
    w::Union{String,UInt16}
    v::Union{String,UInt16}
end
eval(gate::RShiftGate, wires) = get_value(wires, gate.w) >> get_value(wires, gate.v)

function parse_arg(s)
    v = tryparse(Int, s)
    if v === nothing
        return String(s)
    end
    return UInt16(v)
end

function evaluate_circuit(file)
    wires = Dict{String,WireInfo}()
    for line in eachline(file)
        op, name = split(line, " -> ")
        m = match(r"NOT ([a-z]+|\d+)", op)
        if m !== nothing
            wires[name] = WireInfo(NotGate(parse_arg(m[1])))
            continue
        end
        m = match(r"([a-z]+|\d+) AND ([a-z]+|\d+)", op)
        if m !== nothing
            wires[name] = WireInfo(AndGate(parse_arg(m[1]), parse_arg(m[2])))
            continue
        end
        m = match(r"([a-z]+|\d+) OR ([a-z]+|\d+)", op)
        if m !== nothing
            wires[name] = WireInfo(OrGate(parse_arg(m[1]), parse_arg(m[2])))
            continue
        end
        m = match(r"([a-z]+|\d+) LSHIFT ([a-z]+|\d+)", op)
        if m !== nothing
            wires[name] = WireInfo(LShiftGate(parse_arg(m[1]), parse_arg(m[2])))
            continue
        end
        m = match(r"([a-z]+|\d+) RSHIFT ([a-z]+|\d+)", op)
        if m !== nothing
            wires[name] = WireInfo(RShiftGate(parse_arg(m[1]), parse_arg(m[2])))
            continue
        end
        m = match(r"^([a-z]+|\d+)$", op)
        if m !== nothing
            wires[name] = WireInfo(Value(parse_arg(m[1])))
            continue
        end
        error(op)
    end
    return Int(get_value(wires, "a"))
end

@show evaluate_circuit(ARGS[1])

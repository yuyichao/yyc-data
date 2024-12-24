#!/usr/bin/julia

struct Gate
    ty::Symbol
    output::String
    in1::String
    in2::String
end

function apply_renames(gates, renames)
    for i in 1:length(gates)
        gate = gates[i]
        in1 = get(renames, gate.in1, gate.in1)
        in2 = get(renames, gate.in2, gate.in2)
        output = get(renames, gate.output, gate.output)
        if in1 < in2
            gates[i] = Gate(gate.ty, output, in1, in2)
        else
            gates[i] = Gate(gate.ty, output, in2, in1)
        end
    end
end

function rename1(gates)
    renames = Dict{String,String}()
    for gate in gates
        m = match(r"x([0-9][0-9])-y\1", gate.in1 * "-" * gate.in2)
        if m === nothing
            continue
        end
        if gate.output[1] == 'z'
            continue
        end
        if gate.ty == :AND
            renames[gate.output] = "H" * m[1]
        elseif gate.ty == :XOR
            renames[gate.output] = "B" * m[1]
        else
            error("...")
        end
    end
    apply_renames(gates, renames)
end

function rename2(gates)
    renames = Dict{String,String}()
    for gate in gates
        m = match(r"B([0-9][0-9])-z\1", gate.in1 * "-" * gate.output)
        if m === nothing
            continue
        end
        if gate.in2[1] == 'z'
            continue
        end
        if gate.ty == :XOR
            renames[gate.in2] = "C" * m[1]
        end
    end
    apply_renames(gates, renames)
end

function rename3(gates)
    renames = Dict{String,String}()
    for gate in gates
        m = match(r"B([0-9][0-9])-C\1", gate.in1 * "-" * gate.in2)
        if m === nothing
            continue
        end
        if gate.output[1] == 'z'
            continue
        end
        if gate.ty == :AND
            renames[gate.output] = "D" * m[1]
        else
            error("...")
        end
    end
    apply_renames(gates, renames)
end

function eval_gate(file)
    gates = Gate[]
    for line in eachline(file)
        m = match(r"([^ ]*) (AND|OR|XOR) ([^ ]*) -> ([^ ]*)", line)
        if m[1] < m[3]
            push!(gates, Gate(Symbol(m[2]), m[4], m[1], m[3]))
        else
            push!(gates, Gate(Symbol(m[2]), m[4], m[3], m[1]))
        end
    end
    rename1(gates)
    rename2(gates)
    rename3(gates)
    for gate in gates
        println("$(gate.in1) $(gate.ty) $(gate.in2) -> $(gate.output)")
    end
end

eval_gate(ARGS[1])

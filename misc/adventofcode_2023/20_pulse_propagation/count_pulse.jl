#!/usr/bin/julia

mutable struct FlipFlop
    on::Bool
    const outputs::Vector{String}
    FlipFlop(outputs) = new(false, outputs)
end

function add_input!(::FlipFlop, name)
end

function input_pulse(mod::FlipFlop, mod_name, input_name, value, modules)
    if value
        return
    end
    mod.on = !mod.on
    queue_pulse(modules, mod_name, mod.outputs, mod.on)
    return
end

struct Conjunction
    inputs::Dict{String,Bool}
    outputs::Vector{String}
    Conjunction(outputs) = new(Dict{String,Bool}(), outputs)
end

function add_input!(mod::Conjunction, name)
    mod.inputs[name] = false
    return
end

function input_pulse(mod::Conjunction, mod_name, input_name, value, modules)
    mod.inputs[input_name] = value
    return queue_pulse(modules, mod_name, mod.outputs, !all(values(mod.inputs)))
end

const PulseModule = Union{FlipFlop,Conjunction}

struct Pulse
    from::String
    to::String
    value::Bool
end

struct Modules
    broadcast::Vector{String}
    modules::Dict{String,PulseModule}
    pulses::Vector{Pulse}
    Modules() = new(String[], Dict{String,PulseModule}(), Pulse[])
end

function queue_pulse(modules::Modules, mod_name, outputs, value)
    for output in outputs
        push!(modules.pulses, Pulse(mod_name, output, value))
    end
end

function process_pulses(modules::Modules)
    lo = 0
    hi = 0
    while !isempty(modules.pulses)
        pulse = popfirst!(modules.pulses)
        if pulse.value
            hi += 1
        else
            lo += 1
        end
        if haskey(modules.modules, pulse.to)
            input_pulse(modules.modules[pulse.to], pulse.to,
                        pulse.from, pulse.value, modules)
        end
    end
    return (lo, hi)
end

function process_button(modules::Modules)
    queue_pulse(modules, "broadcaster", modules.broadcast, false)
    lo, hi = process_pulses(modules)
    return lo + 1, hi
end

function load_modules(file)
    modules = Modules()
    for line in eachline(file)
        m = match(r"([^ ]+) +-> +(.*)", line)
        outputs = strip.(split(m[2], ","))
        mod = m[1]
        if mod[1] == '%'
            mod_name = mod[2:end]
            modules.modules[mod_name] = FlipFlop(outputs)
        elseif mod[1] == '&'
            mod_name = mod[2:end]
            modules.modules[mod_name] = Conjunction(outputs)
        else
            @assert mod == "broadcaster"
            @assert isempty(modules.broadcast)
            append!(modules.broadcast, outputs)
        end
    end
    for name in modules.broadcast
        add_input!(modules.modules[name], "broadcaster")
    end
    for (mod_name, mod) in modules.modules
        for name in mod.outputs
            if haskey(modules.modules, name)
                add_input!(modules.modules[name], mod_name)
            end
        end
    end
    return modules
end

function count_pulse(file)
    modules = load_modules(file)
    counts = (0, 0)
    for i in 1:1000
        counts = counts .+ process_button(modules)
    end
    return counts
end

@show count_pulse(ARGS[1])

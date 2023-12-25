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

get_state(mod::FlipFlop) = mod.on

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

get_state(mod::Conjunction) = !all(values(mod.inputs))

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

function process_pulses(modules::Modules, checkmod)
    nhi = 0
    while !isempty(modules.pulses)
        pulse = popfirst!(modules.pulses)
        if pulse.value && pulse.from == checkmod
            nhi += 1
        end
        if haskey(modules.modules, pulse.to)
            input_pulse(modules.modules[pulse.to], pulse.to,
                        pulse.from, pulse.value, modules)
        end
    end
    return nhi
end

function process_button(modules::Modules, checkmod)
    queue_pulse(modules, "broadcaster", modules.broadcast, false)
    return process_pulses(modules, checkmod)
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

function collect_gates(modules, startmod, endmod, modset=Set{String}())
    if startmod == endmod
        return modset
    end
    if startmod in modset
        return modset
    end
    push!(modset, startmod)
    for mod in modules.modules[startmod].outputs
        collect_gates(modules, mod, endmod, modset)
    end
    return modset
end

function modset_info(modules, modset, starts, ends)
    startmod = first(intersect(modset, starts))
    endmod = first(intersect(modset, ends))

    modset = sort!(collect(modset))
    modset_states = Dict{Vector{Bool},Int}()

    for mod in modset
        m = modules.modules[mod]
        if isa(m, FlipFlop)
            m.on = false
        end
    end

    step = 0
    while true
        states = [get_state(modules.modules[mod]) for mod in modset]
        if haskey(modset_states, states)
            return step, modset_states[states]
        end
        modset_states[states] = step
        queue_pulse(modules, "broadcaster", [startmod], false)
        np = process_pulses(modules, endmod)
        step += 1
        if np != 0
            @show np, step, get_state(modules.modules[endmod])
        end
    end
end

function count_pulse(file)
    modules = load_modules(file)
    modsets = [collect_gates(modules, mod, "sq") for mod in modules.broadcast]
    nmodsets = length(modsets)
    for i in 2:nmodsets
        for j in 1:i - 1
            @assert isempty(intersect(modsets[i], modsets[j]))
        end
    end

    for modset in modsets
        @show modset_info(modules, modset,
                          modules.broadcast, keys(modules.modules["sq"].inputs))
    end

    # counts = (0, 0)
    # for i in 1:1000
    #     counts = counts .+ process_button(modules)
    # end
    # return counts
end

@show count_pulse(ARGS[1])

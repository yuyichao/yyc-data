#!/usr/bin/julia

using DataStructures

mutable struct Step
    name::String
    ndeps::Int
    dependents::Vector{Step}
end

function sort_steps(file)
    steps = SortedDict{String,Step}()
    get_step(name) = get!(()->Step(name, 0, Step[]), steps, name)

    for line in eachline(file)
        m = match(r"Step (.) must be finished before step (.) can begin\.", line)
        step_before = get_step(m[1])
        step_after = get_step(m[2])
        push!(step_before.dependents, step_after)
        step_after.ndeps += 1
    end

    ready = SortedDict{String,Step}()
    for (name, step) in steps
        if step.ndeps == 0
            ready[name] = step
            delete!(steps, name)
        end
    end

    order = Char[]

    while !isempty(ready)
        name, step = first(ready)
        delete!(ready, name)
        push!(order, name[1])
        for dependent in step.dependents
            dependent.ndeps -= 1
            if dependent.ndeps == 0
                ready[dependent.name] = dependent
                delete!(steps, dependent.name)
            end
        end
    end
    return String(order)
end

@show sort_steps(ARGS[1])

#!/usr/bin/julia

using DataStructures

mutable struct Step
    name::Char
    ndeps::Int
    dependents::Vector{Step}
end

function count_time(file)
    steps = SortedDict{Char,Step}()
    get_step(name) = get!(()->Step(name, 0, Step[]), steps, name)

    for line in eachline(file)
        m = match(r"Step (.) must be finished before step (.) can begin\.", line)
        step_before = get_step(m[1][1])
        step_after = get_step(m[2][1])
        push!(step_before.dependents, step_after)
        step_after.ndeps += 1
    end

    ready = SortedDict{Char,Step}()
    for (name, step) in steps
        if step.ndeps == 0
            ready[name] = step
            delete!(steps, name)
        end
    end

    finish_times = SortedDict{Int,Vector{Step}}()
    max_time = 0

    while !isempty(ready) && length(finish_times) < 5
        name, step = first(ready)
        delete!(ready, name)
        ftime = 60 + (name - 'A') + 1
        finish_times[ftime] = [step]
        @assert max_time < ftime
        max_time = ftime
    end

    while !isempty(finish_times)
        ftime, fsteps = first(finish_times)
        delete!(finish_times, ftime)
        for step in fsteps
            for dependent in step.dependents
                dependent.ndeps -= 1
                if dependent.ndeps == 0
                    ready[dependent.name] = dependent
                    delete!(steps, dependent.name)
                end
            end
        end
        cur_time = ftime
        while !isempty(ready) && length(finish_times) < 5
            name, step = first(ready)
            delete!(ready, name)
            ftime = cur_time + 60 + (name - 'A') + 1
            push!(get!(Vector{Step}, finish_times, ftime), step)
            max_time = max(max_time, ftime)
        end
    end
    return max_time
end

@show count_time(ARGS[1])

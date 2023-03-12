#!/usr/bin/julia

include("utils.jl")

const infile = ARGS[1]

const doc = readxml(infile)
const root = doc.root

const measures = MeasureInfo.(findall("part/measure", root))

mutable struct CostCounter
    single::Int
    double::Int
    group::Int
    fine::Int
    failed::Int
    max_split::Int
    function CostCounter()
        return new(0, 0, 0, 0, 0, 1)
    end
end

function map_measure(info, mid, counter, offset)
    pitches = info.pitches .+ offset
    single = try_map_single(pitches)
    if single === nothing
        double = try_map_double(pitches, info.times, info.total_time, info.group_start)
        if double === nothing
            group = try_map_group(pitches, info.group_start)
            if group === nothing
                fine = try_map_group_fine(pitches, info.times)
                if fine === nothing
                    counter.failed += 1
                else
                    counter.max_split = max(counter.max_split, length(fine))
                    counter.fine += 1
                end
            else
                counter.max_split = max(counter.max_split, length(group))
                counter.group += 1
            end
        else
            counter.max_split = max(counter.max_split, 2)
            counter.double += 1
        end
    else
        counter.single += 1
    end
end

function try_map_all(measures, offset=0)
    counter = CostCounter()
    for (mid, measure) in enumerate(measures)
        map_measure(measure, mid, counter, offset)
    end
    @show counter
end

const all_pitches = collect_all_pitches(root)

const min_offset = min_avail_pitch - all_pitches[1]
const max_offset = max_avail_pitch - all_pitches[end]

if min_offset > max_offset
    error("Pitch range too wide")
end

@time for offset in min_offset:max_offset
    @show offset
    try_map_all(measures, offset)
end

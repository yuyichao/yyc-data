#!/usr/bin/julia

include("utils.jl")

const infile = ARGS[1]

const doc = readxml(infile)
const root = doc.root

const measures = findall("part/measure", root)

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

function map_measure(measure, mid, counter, offset)
    notes = findall("note", measure)
    pitches = Int[]
    times = Int[]
    group_start = Int[]
    cur_t = 0
    for note in notes
        pitch = findfirst("pitch", note)
        if pitch !== nothing
            push!(pitches, parse_pitch_node(pitch) + offset)
        end
        idx = length(pitches)
        is_chord = false
        if findfirst("chord", note) !== nothing
            is_chord = true
            if pitch !== nothing
                push!(times, times[end])
            end
        else
            if pitch !== nothing
                push!(times, cur_t)
            end
            cur_t += parse(Int, nodecontent(findfirst("duration", note)))
        end
        if pitch !== nothing
            beam = findfirst("beam", note)
            tied = findfirst("notations/tied", note)
            if (beam === nothing || strip(nodecontent(beam)) == "begin") &&
                (tied === nothing || tied["type"] == "start")
                push!(group_start, idx)
            end
        end
    end
    single = try_map_single(pitches)
    if single === nothing
        double = try_map_double(pitches, times, cur_t, group_start)
        if double === nothing
            group = try_map_group(pitches, group_start)
            if group === nothing
                fine = try_map_group_fine(pitches, times)
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

for offset in min_offset:max_offset
    @show offset
    try_map_all(measures, offset)
end

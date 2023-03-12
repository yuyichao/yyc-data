#!/usr/bin/julia

include("utils.jl")

const infile = ARGS[1]

const doc = readxml(infile)
const root = doc.root

const measures = findall("part/measure", root)

function collect_all_pitches(root)
    pitches = Set{Int}()
    for pitch in findall("//pitch", root)
        push!(pitches, parse_pitch_node(pitch))
    end
    return sort!(collect(pitches))
end

const all_pitches = collect_all_pitches(root)

# function cost_func(option)
#     min_press = 7
#     for opt in option
#         str_id, pos = opt
#         if pos == 0
#             continue
#         end
#         if str_id < min_press
#             min_press = str_id
#         end
#     end
#     min_use = 7
#     for opt in option
#         str_id, pos = opt
#         if str_id < min_use
#             min_use = str_id
#         end
#     end
#     return (min_press, min_use)
# end

function try_map_all(offset=0)
    single_count = 0
    double_count = 0
    group_count = 0
    failed_count = 0
    for (mid, measure) in enumerate(measures)
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
        # @show group_start
        # @show cur_t
        single = try_map_stable(pitches)
        # @show single
        if single === nothing
            double = try_map_stable2_mid(pitches, times, cur_t, group_start)
            if double === nothing
                group = try_map_stable_n(pitches, group_start)
                if group === nothing
                    failed_count += 1
                else
                    group_count += 1
                end
            else
                double_count += 1
            end
        else
            single_count += 1
        end
    end
    @show single_count, double_count, group_count, failed_count
end

try_map_all(-3)

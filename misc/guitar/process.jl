#!/usr/bin/julia

using EzXML

const infile = ARGS[1]

const doc = readxml(infile)
const root = doc.root

function pitch_to_num(step, octave, alter=0)
    if step == "C"
        num = 0
    elseif step == "D"
        num = 2
    elseif step == "E"
        num = 4
    elseif step == "F"
        num = 5
    elseif step == "G"
        num = 7
    elseif step == "A"
        num = 9
    elseif step == "B"
        num = 11
    else
        error("Invalid step: $step")
    end
    num += octave * 12
    return num + alter
end

function num_to_pitch(num)
    octave = num ÷ 12
    num = num % 12
    if num == 0
        return ("C", octave, 0)
    elseif num == 1
        return ("C", octave, 1)
    elseif num == 2
        return ("D", octave, 0)
    elseif num == 3
        return ("D", octave, 1)
    elseif num == 4
        return ("E", octave, 0)
    elseif num == 5
        return ("F", octave, 0)
    elseif num == 6
        return ("F", octave, 1)
    elseif num == 7
        return ("G", octave, 0)
    elseif num == 8
        return ("G", octave, 1)
    elseif num == 9
        return ("A", octave, 0)
    elseif num == 10
        return ("A", octave, 1)
    elseif num == 11
        return ("B", octave, 0)
    end
end

const measures = findall("part/measure", root)

function collect_all_pitches(root)
    pitches = Set{Int}()
    for pitch in findall("//pitch", root)
        step = strip(nodecontent(findfirst("step", pitch)))
        octave = parse(Int, nodecontent(findfirst("octave", pitch)))
        alter = findfirst("alter", pitch)
        alter = alter === nothing ? 0 : parse(Int, nodecontent(alter))
        num = pitch_to_num(step, octave, alter)
        push!(pitches, num)
    end
    return sort!(collect(pitches))
end

const all_pitches = collect_all_pitches(root)

const base_pitches = [pitch_to_num("E", 2), pitch_to_num("A", 2),
                      pitch_to_num("D", 3), pitch_to_num("G", 3),
                      pitch_to_num("B", 3), pitch_to_num("E", 4)]
const max_pitch_diff = 21

function _get_all_options()
    options = Vector{NTuple{2,Int}}[]
    for i in 1:(base_pitches[end] + max_pitch_diff)
        opts = NTuple{2,Int}[]
        push!(options, opts)
        for (str_id, base) in enumerate(base_pitches)
            diff = i - base
            if 0 <= diff <= max_pitch_diff
                push!(opts, (str_id, diff))
            end
        end
    end
    return options
end

const pitch_options = _get_all_options()
const _empty_option = NTuple{2,Int}[]

get_pitch_options(pitch) = get(pitch_options, pitch, _empty_option)

# for measure in measures
#     notes = findall("note", measure)
#     for note in notes
#         pitch = findfirst("pitch", note)
#         println(pitch)
#     end
# end

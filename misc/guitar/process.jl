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

function parse_pitch_node(pitch)
    step = strip(nodecontent(findfirst("step", pitch)))
    octave = parse(Int, nodecontent(findfirst("octave", pitch)))
    alter = findfirst("alter", pitch)
    alter = alter === nothing ? 0 : parse(Int, nodecontent(alter))
    return pitch_to_num(step, octave, alter)
end

function collect_all_pitches(root)
    pitches = Set{Int}()
    for pitch in findall("//pitch", root)
        push!(pitches, parse_pitch_node(pitch))
    end
    return sort!(collect(pitches))
end

const all_pitches = collect_all_pitches(root)

const base_pitches = [pitch_to_num("E", 2), pitch_to_num("A", 2),
                      pitch_to_num("D", 3), pitch_to_num("G", 3),
                      pitch_to_num("B", 3), pitch_to_num("E", 4)]
const max_pitch_diff = 21
const max_avail_pitch = base_pitches[end] + max_pitch_diff

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

function find_cover_nth!(all_res, current, used, options, n)
    if n > length(options)
        push!(all_res, copy(current))
        return
    end
    for opt in options[n]
        str_id = opt[1]
        if used[str_id]
            continue
        end
        current[n] = opt
        used[str_id] = true
        find_cover_nth!(all_res, current, used, options, n + 1)
        used[str_id] = false
    end
end

function max_pos_diff(option)
    min_pos = 21
    max_pos = 1
    for opt in option
        pos = opt[2]
        if pos == 0
            continue
        end
        if pos > max_pos
            max_pos = pos
        end
        if pos < min_pos
            min_pos = pos
        end
    end
    if max_pos <= min_pos
        return 0
    end
    return max_pos - min_pos
end

function cost_func(option)
    min_press = 7
    for opt in option
        str_id, pos = opt
        if pos == 0
            continue
        end
        if str_id < min_press
            min_press = str_id
        end
    end
    min_use = 7
    for opt in option
        str_id, pos = opt
        if str_id < min_use
            min_use = str_id
        end
    end
    return (min_press, min_use)
end

# Try to map all the pitches to a stable position that is no more than
# 2-3 positions apart other than 0
function try_map_stable(pitches)
    sorted_pitches = sort!(collect(Set(pitches)))
    options = get_pitch_options.(sorted_pitches)
    all_res = Vector{NTuple{2,Int}}[]
    if any(isempty, options) || length(options) > 6
        return all_res
    end
    nnotes = length(options)
    find_cover_nth!(all_res, Vector{NTuple{2,Int}}(undef, nnotes),
                    falses(6), options, 1)
    all_res = [res for res in all_res if max_pos_diff(res) <= 2]
    if isempty(all_res)
        all_res = [res for res in all_res if max_pos_diff(res) <= 3]
    end
    sort!(all_res, by=cost_func, rev=true)
    idx_map = Dict(pitch=>id for (id, pitch) in enumerate(sorted_pitches))
    return ([[res[idx_map[p]] for p in pitches] for res in all_res],
            [(minimum(r[2] for r in res if r[2] != 0),
              maximum(r[2] for r in res if r[2] != 0)) for res in all_res])
end

function try_map_all(offset=0)
    for measure in measures
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
                if beam === nothing || strip(nodecontent(beam)) == "begin"
                    push!(group_start, idx)
                end
            end
        end
        @show group_start
        @show cur_t
        @show try_map_stable(pitches)
    end
end

try_map_all(-3)

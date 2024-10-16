#!/usr/bin/julia

using EzXML

# Basic conversion functions
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

# Document utilities
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

struct MeasureInfo
    node::EzXML.Node
    notes::Vector{EzXML.Node}
    # The following fields only consider note with non-empty pitch
    # This includes the group start index, which indexes into the pitch
    # and not the notes.
    pitch_noteidx::Vector{Int}
    pitches::Vector{Int}
    times::Vector{Int}
    group_start::Vector{Int}
    total_time::Int
    function MeasureInfo(node)
        notes = findall("note", node)
        pitch_noteidx = Int[]
        pitches = Int[]
        times = Int[]
        group_start = Int[]
        total_time = 0
        for (noteidx, note) in enumerate(notes)
            pitch = findfirst("pitch", note)
            if pitch !== nothing
                push!(pitch_noteidx, noteidx)
                push!(pitches, parse_pitch_node(pitch))
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
                    push!(times, total_time)
                end
                total_time += parse(Int, nodecontent(findfirst("duration", note)))
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
        return new(node, notes, pitch_noteidx, pitches, times,
                   group_start, total_time)
    end
end

# Properties of guitar
const base_pitches = [pitch_to_num("E", 2), pitch_to_num("A", 2),
                      pitch_to_num("D", 3), pitch_to_num("G", 3),
                      pitch_to_num("B", 3), pitch_to_num("E", 4)]
const max_pitch_diff = 21
const max_avail_pitch = base_pitches[end] + max_pitch_diff
const min_avail_pitch = base_pitches[1]

struct NoteMapping
    str::Int
    pos::Int
end

function _get_all_options()
    options = Vector{NoteMapping}[]
    for i in 1:(base_pitches[end] + max_pitch_diff)
        opts = NoteMapping[]
        push!(options, opts)
        for (str_id, base) in enumerate(base_pitches)
            diff = i - base
            if 0 <= diff <= max_pitch_diff
                push!(opts, NoteMapping(str_id, diff))
            end
        end
    end
    return options
end

const pitch_options = _get_all_options()
const _empty_option = NoteMapping[]

get_pitch_options(pitch) = get(pitch_options, pitch, _empty_option)

# For finding note to technique options.
function find_cover_nth!(all_res, current, used, options, n)
    if n > length(options)
        push!(all_res, copy(current))
        return
    end
    for opt in options[n]
        str_id = opt.str
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
        pos = opt.pos
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

function finger_count(option)
    option = sort(option, by=x->(x.pos, x.str))
    count = 0
    used_string = 0
    pressed = 0
    for opt in option
        if opt.str == 0
            used_string = max(used_string, opt.pos)
            continue
        end
        # We are already pressing this position in a row
        if pressed == opt.pos
            used_string = max(used_string, opt.pos)
            continue
        end
        if opt.str > used_string
            pressed = opt.pos
        end
        count += 1
        used_string = max(used_string, opt.pos)
    end
    return count
end

# Try to map all the pitches to a stable position that is no more than
# 2-3 positions apart other than 0
function try_map_single(pitches)
    sorted_pitches = sort!(collect(Set(pitches)))
    options = get_pitch_options.(sorted_pitches)
    all_res = Vector{NoteMapping}[]
    if any(isempty, options) || length(options) > 6
        return
    end
    nnotes = length(options)
    find_cover_nth!(all_res, Vector{NoteMapping}(undef, nnotes),
                    falses(6), options, 1)
    _all_res = [res for res in all_res if max_pos_diff(res) <= 2]
    # if isempty(_all_res)
    #     _all_res = [res for res in all_res if max_pos_diff(res) <= 3]
    # end
    all_res = _all_res
    if isempty(all_res)
        return
    end
    fcs = finger_count.(all_res)
    min_finger = minimum(fcs)
    # all_res = all_res[fcs .<= min_finger + 1]
    all_res = all_res[fcs .<= min_finger]
    idx_map = Dict(pitch=>id for (id, pitch) in enumerate(sorted_pitches))
    return [[res[idx_map[p]] for p in pitches] for res in all_res]
end

function try_map_double(pitches, times, total_t, group_start)
    if length(group_start) < 2
        return
    end
    mid_time = total_t ÷ 2
    min_diff_from_mid = mid_time
    mid_start = 0
    for idx in group_start
        t = times[idx]
        d = abs(t - mid_time)
        if d < min_diff_from_mid
            mid_start = idx
            min_diff_from_mid = d
        end
    end
    res1 = try_map_single(pitches[1:mid_start - 1])
    if res1 !== nothing
        res2 = try_map_single(pitches[mid_start:end])
        if res2 !== nothing
            return (mid_start - 1, length(pitches) - mid_start + 1), (res1, res2)
        end
    end
    return
end

function try_map_group(pitches, group_start)
    lens = Int[]
    res = Vector{Vector{NoteMapping}}[]
    grp_id = 1
    ngrp = length(group_start)
    while grp_id <= ngrp
        found = false
        for i in ngrp:-1:grp_id
            idx1 = group_start[grp_id]
            idx2 = i == ngrp ? length(pitches) : (group_start[i + 1] - 1)
            seg = try_map_single(pitches[idx1:idx2])
            if seg !== nothing
                found = true
                push!(lens, idx2 - idx1 + 1)
                push!(res, seg)
                grp_id = i + 1
                break
            end
        end
        if !found
            return
        end
    end
    return lens, res
end

function try_map_group_fine(pitches, times)
    group_start = Int[1]
    for i in 2:length(pitches)
        if pitches[i] == pitches[i - 1] || times[i] == times[i - 1]
            continue
        end
        push!(group_start, i)
    end
    return try_map_group(pitches, group_start)
end

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

struct NoteID
    mid::Int
    noteidx::Int
    pitchidx::Int
end

struct SegmentMappings
    notes::Vector{NoteID}
    options::Vector{Vector{NoteMapping}}
end

function map_measure(segments, info, mid, counter, offset)
    pitches = info.pitches .+ offset
    noteids = [NoteID(mid, info.pitch_noteidx[i], i) for i in 1:length(pitches)]
    single = try_map_single(pitches)
    if single !== nothing
        counter.single += 1
        push!(segments, SegmentMappings(noteids, single))
        return
    end
    double = try_map_double(pitches, info.times, info.total_time, info.group_start)
    if double !== nothing
        counter.max_split = max(counter.max_split, 2)
        counter.double += 1
        (n1, n2), (res1, res2) = double
        push!(segments, SegmentMappings(noteids[1:n1], res1))
        push!(segments, SegmentMappings(noteids[n1 + 1:end], res2))
        return
    end
    group = try_map_group(pitches, info.group_start)
    if group !== nothing
        lens, group = group
        ngrp = length(group)
        counter.max_split = max(counter.max_split, ngrp)
        counter.group += 1
        istart = 1
        for i in 1:ngrp
            istart_next = istart + lens[i]
            push!(segments, SegmentMappings(noteids[istart:istart_next - 1],
                                            group[i]))
            istart = istart_next
        end
        return
    end
    fine = try_map_group_fine(pitches, info.times)
    if fine !== nothing
        lens, group = fine
        ngrp = length(group)
        counter.max_split = max(counter.max_split, ngrp)
        counter.fine += 1
        istart = 1
        for i in 1:ngrp
            istart_next = istart + lens[i]
            push!(segments, SegmentMappings(noteids[istart:istart_next - 1],
                                            group[i]))
            istart = istart_next
        end
        return
    end
    counter.failed += 1
end

struct MapAllResult
    offset::Int
    counter::CostCounter
    segments::Vector{SegmentMappings}
    MapAllResult(offset) = new(offset, CostCounter(), SegmentMappings[])
end

function map_all_measures(measures, offset=0)
    res = MapAllResult(offset)
    for (mid, measure) in enumerate(measures)
        map_measure(res.segments, measure, mid, res.counter, offset)
    end
    return res
end

#!/usr/bin/julia

include("utils.jl")

const infile = ARGS[1]
const outfile = ARGS[2]

const doc = readxml(infile)
const root = doc.root

const measures = MeasureInfo.(findall("part/measure", root))

const all_pitches = collect_all_pitches(root)

const min_offset = min_avail_pitch - all_pitches[1]
const max_offset = max_avail_pitch - all_pitches[end]

if min_offset > max_offset
    error("Pitch range too wide")
end

const map_all_results = [map_all_measures(measures, offset)
                         for offset in min_offset:max_offset]

sort!(map_all_results, by=x->begin
          c = x.counter
          return (c.failed, c.max_split, -c.single, -c.double, -c.group)
      end)
const map_all_result = map_all_results[1]
@assert map_all_result.counter.failed == 0

function compare_mapping(mapping1, mapping2)
    set1 = Set(mapping1)
    set2 = Set(mapping2)
    overlap = length(union(set1, set2))
    minpos1, maxpos1 = extrema(m.pos for m in mapping1)
    minpos2, maxpos2 = extrema(m.pos for m in mapping2)
    maxdist = max(abs(maxpos2 - minpos1), abs(maxpos1 - minpos2))
    return overlap, maxdist
end

function cost_func(mapping1, mapping2)
    overlap, maxdist = compare_mapping(mapping1, mapping2)
    return maxdist - 2 * overlap
end

function pick_segments(result)
    paths = Tuple{Int,Vector{Int}}[]
    local prev_options
    for seg in result.segments
        options = seg.options
        if isempty(paths)
            # Special case for the first element
            for i in 1:length(options)
                push!(paths, (0, [i]))
            end
            prev_options = options
            continue
        end
        new_paths = Tuple{Int,Vector{Int}}[]
        for (oidx, option) in enumerate(options)
            min_cost = typemax(Int)
            min_pidx = 0
            for (pidx, (prev_cost, path)) in enumerate(paths)
                prev_option = prev_options[pidx]
                cost = prev_cost + cost_func(prev_option, option)
                if cost < min_cost
                    min_cost = cost
                    min_pidx = pidx
                end
            end
            push!(new_paths, (min_cost, [paths[min_pidx][2]; oidx]))
        end
        paths = new_paths
        prev_options = options
    end
    local optimal_path
    if length(paths) <= 1
        optimal_path = paths[1][2]
    else
        min_cost = typemax(Int)
        for (pidx, (cost, path)) in enumerate(paths)
            if cost < min_cost
                min_cost = cost
                optimal_path = path
            end
        end
    end
    for (i, p) in enumerate(optimal_path)
        options = result.segments[i].options
        options[1] = options[p]
        resize!(options, 1)
    end
end

pick_segments(map_all_result)

function set_node_content(parent, name, content)
    node = findfirst(name, parent)
    if node === nothing
        return addelement!(parent, name, content)
    else
        setnodecontent!(node, content)
        return node
    end
end

function update_doc!(measures, result)
    for seg in result.segments
        note_ids = seg.notes
        option = seg.options[1]
        if result.offset != 0
            for (i, (note_id, opt)) in enumerate(zip(note_ids, option))
                measure = measures[note_id.mid]
                origin_pitch = measure.pitches[note_id.pitchidx]
                new_pitch = origin_pitch + result.offset
                step, octave, alter = num_to_pitch(new_pitch)
                note = measure.notes[note_id.noteidx]
                pitch = findfirst("pitch", note)
                @assert pitch !== nothing
                setnodecontent!(findfirst("step", pitch), step)
                setnodecontent!(findfirst("octave", pitch), string(octave))
                alter_node = findfirst("alter", pitch)
                if alter == 0
                    if alter_node !== nothing
                        unlink!(alter_node)
                    end
                else
                    if alter_node === nothing
                        addelement!(pitch, "alter", string(alter))
                    else
                        setnodecontent!(alter_node, string(alter))
                    end
                end
            end
        end
        for (i, (note_id, opt)) in enumerate(zip(note_ids, option))
            measure = measures[note_id.mid]
            note = measure.notes[note_id.noteidx]
            notations = findfirst("notations", note)
            if notations === nothing
                notations = addelement!(note, "notations", "")
            end
            technical = findfirst("technical", notations)
            if technical === nothing
                technical = addelement!(notations, "technical", "")
            end
            set_node_content(technical, "string", string(7 - opt.str))
            set_node_content(technical, "fret", string(opt.pos))
        end
    end
end

update_doc!(measures, map_all_result)

open(outfile, "w") do io
    print(io, doc)
end

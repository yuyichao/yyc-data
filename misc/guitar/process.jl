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

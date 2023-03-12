#!/usr/bin/julia

include("utils.jl")

const infile = ARGS[1]

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

function _pick_segments(result, start_idx, end_idx)
    sz = end_idx - start_idx + 1
    min_idxs = zeros(Int, sz)
    min_cost = typemax(Int)
    all_cost = Vector{Matrix{Int}}(undef, sz + 1)
    lens = [length(result.segments[i].options) for i in start_idx:end_idx]
    for i in 1:sz + 1
        sz1 = i == 1 ? 1 : lens[i - 1]
        sz2 = i == sz + 1 ? 1 : lens[i]
        i1 = i + start_idx - 2
        i2 = i + start_idx - 1
        if i1 < 1 || i2 > length(result.segments)
            all_cost[i] = zeros(Int, sz1, sz2)
            continue
        end
        cost_mat = Matrix{Int}(undef, sz1, sz2)
        all_cost[i] = cost_mat
        for idx2 in 1:sz2
            for idx1 in 1:sz1
                cost_mat[idx1, idx2] = cost_func(result.segments[i1].options[idx1],
                                                 result.segments[i2].options[idx2])
            end
        end
    end
    if true
        # Brute force search...
        # Good enough for now...
        idxs = ones(Int, sz)
        while all(idxs .<= lens)
            cost = 0
            for i in 1:sz + 1
                idx1 = i == 1 ? 1 : idxs[i - 1]
                idx2 = i == sz + 1 ? 1 : idxs[i]
                cost += all_cost[i][idx1, idx2]
            end
            if cost < min_cost
                min_idxs .= idxs
                min_cost = cost
            end
            for i in sz:-1:1
                idxs[i] += 1
                if idxs[i] <= lens[i]
                    break
                end
                if i == 1
                    break
                end
                idxs[i] = 1
            end
        end
        for i in 1:sz
            idx = min_idxs[i]
            options = result.segments[i + start_idx - 1].options
            options[1] = options[idx]
            resize!(options, 1)
        end
    end
end

function pick_segments(result)
    start_idx = 0
    for (i, seg) in enumerate(result.segments)
        if length(seg.options) == 1
            if start_idx == 0
                continue
            end
            _pick_segments(result, start_idx, i - 1)
            start_idx = 0
            continue
        end
        if start_idx == 0
            start_idx = i
        end
    end
    if start_idx != 0
        _pick_segments(result, start_idx, length(result.segments))
        start_idx = 0
    end
end

pick_segments(map_all_result)

@show [length(segment.options) for segment in map_all_result.segments]
@show prod(big(length(segment.options)) for segment in map_all_result.segments)

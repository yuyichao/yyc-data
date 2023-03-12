#!/usr/bin/julia

include("utils.jl")

const infile = ARGS[1]

const doc = readxml(infile)
const root = doc.root

const measures = MeasureInfo.(findall("part/measure", root))

function try_map_all(measures, offset=0)
    counter = CostCounter()
    segments = SegmentMappings[]
    for (mid, measure) in enumerate(measures)
        map_measure(segments, measure, mid, counter, offset)
    end
    # @show counter
end

const all_pitches = collect_all_pitches(root)

const min_offset = min_avail_pitch - all_pitches[1]
const max_offset = max_avail_pitch - all_pitches[end]

if min_offset > max_offset
    error("Pitch range too wide")
end

@time for offset in min_offset:max_offset
    # @show offset
    try_map_all(measures, offset)
end

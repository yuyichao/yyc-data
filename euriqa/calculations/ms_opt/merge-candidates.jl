#!/usr/bin/julia

include("opt-all-pair.jl")

const prefix = ARGS[1]
const block_size = 2000
meta, candidates = load_candidates_dir(dirname(prefix))
tail_candidates = Candidate[]
function cut_tail!(candidates, tail_candidates)
    n = length(candidates)
    nround = n ÷ block_size * block_size
    append!(tail_candidates, @view(candidates[nround + 1:end]))
    resize!(candidates, nround)
    return
end
cut_tail!(candidates, tail_candidates)
println("Load original: $(length(candidates)) [tail: $(length(tail_candidates))]")
for d in ARGS[2:end]
    load_candidates_dir(d, meta=meta, candidates=candidates)
    cut_tail!(candidates, tail_candidates)
    println("Loaded $d: $(length(candidates)) [tail: $(length(tail_candidates))]")
end
append!(candidates, tail_candidates)
println("Total $(length(candidates))")
save_candidates(prefix, candidates, meta, block_size=block_size, use_pb=true)

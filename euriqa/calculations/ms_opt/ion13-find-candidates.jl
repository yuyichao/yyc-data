#!/usr/bin/julia

include("opt-all-pair.jl")
include("ion13-params.jl")

const prefix = ARGS[1]
const nrep = parse(Int, ARGS[2])

meta = Dict("amp_ratio"=>amp_ratio, "nseg"=>nseg)
_meta, candidates = load_candidates_dir(dirname(prefix))
println("Loaded $(length(candidates))")
@assert _meta === nothing || _meta == meta
const pre_pool = ThreadObjectPool() do
    return PreOptimizer{nseg}(2π .* fs, 2π .* (fs .+ 0.3); amp_ratio=amp_ratio,
                              tmin=350, tmax=600, ntimes=21,
                              ωmin=2π * 2.28, ωmax=2π * 2.497)
end
candidates = @time opt_all_rounds!(pre_pool, nrep, candidates)
@show length(candidates)
save_candidates(prefix, candidates, meta)

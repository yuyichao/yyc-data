#!/usr/bin/julia

include("multi-pairs-test-utils.jl")

const nsegs = [80, 100, 120, 140, 160, 180, 200]
const compute_bufs_opt = [SL.ComputeBuffer{nseg,Float64}(Val(SS.mask_allδ), Val(SS.mask_allδ)) for nseg in nsegs]
const preobjs = [get_preobj(modes, buf, 0.6) for buf in compute_bufs_opt]
const preopts = [PreOptimizer(preobj) for preobj in preobjs]

for i in 1:length(nsegs)
    precompile(run_preopts, (typeof(nsegs[i]), typeof(preopts[i]),
                             Int, Float64, Float64, Int))
end

@time Threads.@threads for i in 1:length(nsegs)
    open("pairs_data_$(nsegs[i])_segs_3.binpb", "w") do io
        encoder = PB.ProtoEncoder(io)
        PB.encode(encoder, run_preopts(nsegs[i], preopts[i], 50, 150.0, 350.0, 25))
    end
end

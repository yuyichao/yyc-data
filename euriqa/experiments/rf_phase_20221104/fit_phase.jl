#!/usr/bin/julia

include("fitting_utils.jl")

const data_1 = load_data(joinpath(@__DIR__, "data/C2Trace00000.h5"))
const blocks_1 = find_blocks(data_1)
const tvblocks_1 = get_tvblocks(data_1, blocks_1)

uncs1 = fit_multi_phases(tvblocks_1...).uncs[4:end]
uncs2 = fit_multi_phases(shrink_block.(tvblocks_1[1], 0.1, 0.9),
                         shrink_block.(tvblocks_1[2], 0.1, 0.9)).uncs[4:end]
uncs3 = fit_multi_phases(shrink_block.(tvblocks_1[1], 0.3, 0.7),
                         shrink_block.(tvblocks_1[2], 0.3, 0.7)).uncs[4:end]

# @show shrink_block.(find_blocks(load_data(ARGS[1])), 0.1, 0.9)

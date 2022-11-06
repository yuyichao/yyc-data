#!/usr/bin/julia

include("fitting_utils.jl")

const names_C2 = ["C2Trace00000.h5", "C2Trace00001.h5", "C2Trace00002.h5",
                  "C2Trace00003.h5", "C2Trace00005.h5", "C2Trace00006.h5",
                  "C2Trace00007.h5", "C2Trace00008.h5", "C2Trace00009.h5",
                  "C2Trace00010.h5"]
const names_C3 = ["C3Trace00000.h5", "C3Trace00001.h5", "C3Trace00002.h5",
                  "C3Trace00003.h5", "C3Trace00004.h5", "C3Trace00005.h5",
                  "C3Trace00007.h5", "C3Trace00008.h5", "C3Trace00009.h5",
                  "C3Trace00010.h5"]

const data_set_C2 = [load_data(joinpath(@__DIR__, "data", name)) for name in names_C2]
const data_set_C3 = [load_data(joinpath(@__DIR__, "data", name)) for name in names_C3]

struct BlockInfo
    data_idx::Int
    block_idx::Int
    time_offset::Float64
end

function collect_blocks(data_set, r_start=0, r_end=1)
    block_infos = BlockInfo[]
    tblocks = Vector{Float64}[]
    vblocks = Vector{Float64}[]
    for (data_idx, data) in enumerate(data_set)
        idx_blocks = find_blocks(data)
        for (block_idx, (tb, vb)) in enumerate(zip(get_tvblocks(data, idx_blocks)...))
            toffset = (tb[1] + tb[end]) / 2
            push!(tblocks, shrink_block(tb .- toffset, r_start, r_end))
            push!(vblocks, shrink_block(vb, r_start, r_end))
            push!(block_infos, BlockInfo(data_idx, block_idx, toffset))
        end
    end
    return block_infos, tblocks, vblocks
end

const block_infos_C2, tblocks_C2, vblocks_C2 = collect_blocks(data_set_C2)
const block_infos_C3, tblocks_C3, vblocks_C3 = collect_blocks(data_set_C3)
const fit_C2 = fit_multi_phases(tblocks_C2, vblocks_C2)
const fit_C3 = fit_multi_phases(tblocks_C3, vblocks_C3)

# const block_infos_sC2, tblocks_sC2, vblocks_sC2 = collect_blocks(data_set_C2, 0.2, 0.8)
# const block_infos_sC3, tblocks_sC3, vblocks_sC3 = collect_blocks(data_set_C3, 0.2, 0.8)
# const fit_sC2 = fit_multi_phases(tblocks_sC2, vblocks_sC2)
# const fit_sC3 = fit_multi_phases(tblocks_sC3, vblocks_sC3)

# @show fit_sC2.uncs .- fit_C2.uncs
# @show fit_sC3.uncs .- fit_C3.uncs

@show fit_C2.uncs[3]
@show fit_C2.uncs[4:end]
@show fit_C3.uncs[3]
@show fit_C3.uncs[4:end]

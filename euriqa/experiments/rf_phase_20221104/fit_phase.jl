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
    block_counts = Int[]
    block_infos = BlockInfo[]
    tblocks = Vector{Float64}[]
    vblocks = Vector{Float64}[]
    for (data_idx, data) in enumerate(data_set)
        idx_blocks = find_blocks(data)
        push!(block_counts, length(idx_blocks))
        for (block_idx, (tb, vb)) in enumerate(zip(get_tvblocks(data, idx_blocks)...))
            toffset = (tb[1] + tb[end]) / 2
            push!(tblocks, shrink_block(tb .- toffset, r_start, r_end))
            push!(vblocks, shrink_block(vb, r_start, r_end))
            push!(block_infos, BlockInfo(data_idx, block_idx, toffset))
        end
    end
    return block_infos, block_counts, tblocks, vblocks
end

const block_infos_C2, block_counts_C2, tblocks_C2, vblocks_C2 =
    collect_blocks(data_set_C2)
const block_infos_C3, block_counts_C3, tblocks_C3, vblocks_C3 =
    collect_blocks(data_set_C3)
const fit_C2 = fit_multi_phases(tblocks_C2, vblocks_C2)
const fit_C3 = fit_multi_phases(tblocks_C3, vblocks_C3)

function find_block_info(block_infos, data_idx, block_idx)
    function lt(block_info::BlockInfo, idx::Tuple)
        return (block_info.data_idx, block_info.block_idx) < idx
    end
    function lt(idx::Tuple, block_info::BlockInfo)
        return idx < (block_info.data_idx, block_info.block_idx)
    end
    rng = searchsorted(block_infos, (data_idx, block_idx), lt=lt)
    if length(rng) == 0
        error("Unable to find block info")
    elseif length(rng) != 1
        error("Duplicated find block info")
    end
    idx = first(rng)
    return idx, block_infos[idx]
end

function diff_phase(block_infos, fit, data_idx, block_idx1, block_idx2)
    idx1, block_info1 = find_block_info(block_infos, data_idx, block_idx1)
    idx2, block_info2 = find_block_info(block_infos, data_idx, block_idx2)

    # The phase (in unit of 2pi) for the two blocks are
    # `f * (t - toffset_1) + ϕf_1` and `f * (t - toffset_2) + ϕf_2`
    # Phase difference is `f * (t - toffset_2) + ϕf_2 - f * (t - toffset_1) - ϕf_1`
    # or `f * (toffset_1 - toffset_2) + ϕf_2 - ϕf_1`
    Δϕf = mod(fit.param[3] * (block_info1.time_offset - block_info2.time_offset), 1)
    Δϕf = mod(Δϕf + fit.param[3 + idx2] - fit.param[3 + idx1] + 0.5, 1) - 0.5
    vec = zeros(length(fit.param))
    vec[3] = block_info1.time_offset - block_info2.time_offset
    vec[3 + idx2] += 1
    vec[3 + idx1] += -1

    s_Δϕf = vec' * fit.covar * vec
    return Unc(Δϕf, sqrt(s_Δϕf))
end

# const block_infos_sC2, block_counts_sC2, tblocks_sC2, vblocks_sC2 =
#     collect_blocks(data_set_C2, 0.1, 0.9)
# const block_infos_sC3, block_counts_sC3, tblocks_sC3, vblocks_sC3 =
#     collect_blocks(data_set_C3, 0.1, 0.9)
# const fit_sC2 = fit_multi_phases(tblocks_sC2, vblocks_sC2)
# const fit_sC3 = fit_multi_phases(tblocks_sC3, vblocks_sC3)
# # @show fit_sC2.uncs .- fit_C2.uncs
# # @show fit_sC3.uncs .- fit_C3.uncs

# function compare_fits(ndata, block_counts, block_infos1, fit1, block_infos2, fit2)
#     for data_idx in 1:ndata
#         nblocks = block_counts[data_idx]
#         for blk_id2 in 1:nblocks
#             for blk_id1 in (blk_id2 + 1):nblocks
#                 d1 = diff_phase(block_infos1, fit1, data_idx, blk_id1, blk_id2)
#                 d2 = diff_phase(block_infos2, fit2, data_idx, blk_id1, blk_id2)
#                 diff = d1 - d2
#                 @show (mod(diff.a + 0.5, 1) - 0.5) / diff.s
#             end
#         end
#     end
# end
# compare_fits(length(names_C2), block_counts_C2, block_infos_C2, fit_C2,
#              block_infos_sC2, fit_sC2)
# compare_fits(length(names_C3), block_counts_C3, block_infos_C3, fit_C3,
#              block_infos_sC3, fit_sC3)

const prefix = joinpath(@__DIR__, "data/")

function save_tables(prefix, names, block_infos, block_counts, fit)
    for (data_idx, name) in enumerate(names)
        name = replace(name, r"\.[^.]*$"=>"")
        open("$(prefix)_$(name).csv", "w") do io
            nblocks = block_counts[data_idx]
            for blk_id1 in 1:nblocks
                print(io, ",$blk_id1")
            end
            println(io)
            for blk_id2 in 1:nblocks
                print(io, blk_id2)
                for blk_id1 in 1:nblocks
                    if blk_id1 <= blk_id2
                        print(io, ",-")
                        continue
                    end
                    d = diff_phase(block_infos, fit, data_idx, blk_id1, blk_id2)
                    print(io, ",$d")
                end
                println(io)
            end
        end
    end
end

save_tables("$(prefix)phase_diff", names_C2, block_infos_C2, block_counts_C2, fit_C2)
save_tables("$(prefix)phase_diff", names_C3, block_infos_C3, block_counts_C3, fit_C3)

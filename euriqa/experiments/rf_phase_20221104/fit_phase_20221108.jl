#!/usr/bin/julia

include("fitting_utils.jl")

using PyPlot
using NaCsPlot

const names_C1 = ["C1Trace00000.h5", "C1Trace00001.h5", "C1Trace00002.h5",
                  "C1Trace00003.h5", "C1Trace00004.h5", "C1Trace00005.h5",
                  "C1Trace00006.h5", "C1Trace00007.h5", "C1Trace00008.h5",
                  "C1Trace00009.h5"]
const names_C2 = ["C2Trace00011.h5", "C2Trace00012.h5", "C2Trace00013.h5",
                  "C2Trace00014.h5", "C2Trace00015.h5", "C2Trace00014.h5",
                  "C2Trace00017.h5", "C2Trace00018.h5", "C2Trace00019.h5",
                  "C2Trace00020.h5"]

const data_set_C1 = [load_data(joinpath(@__DIR__, "data", name)) for name in names_C1]
const data_set_C2 = [load_data(joinpath(@__DIR__, "data", name)) for name in names_C2]

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
            # toffset = (tb[1] + tb[end]) / 2
            # This time offset produce a more linear-phase-shift-free fit.
            # There isn't anything physical about such a shift and I'm too lazy
            # to figure out a way get rid of it more reliably...
            toffset = tb[1]
            push!(tblocks, shrink_block(tb .- toffset, r_start, r_end))
            push!(vblocks, shrink_block(vb, r_start, r_end))
            push!(block_infos, BlockInfo(data_idx, block_idx, toffset))
        end
    end
    return block_infos, block_counts, tblocks, vblocks
end

const block_infos_C1, block_counts_C1, tblocks_C1, vblocks_C1 =
    collect_blocks(data_set_C1)
const block_infos_C2, block_counts_C2, tblocks_C2, vblocks_C2 =
    collect_blocks(data_set_C2)
const fit_C1 = fit_multi_phases(tblocks_C1, vblocks_C1)
const fit_C2 = fit_multi_phases(tblocks_C2, vblocks_C2)

function find_block_info(block_infos, data_idx, block_idx)
    function lt(block_info::BlockInfo, idx::Tuple)
        return (block_info.data_idx, block_info.block_idx) < idx
    end
    function lt(idx::Tuple, block_info::BlockInfo)
        return idx < (block_info.data_idx, block_info.block_idx)
    end
    rng = searchsorted(block_infos, (data_idx, block_idx), lt=lt)
    if length(rng) == 0
        @show data_idx, block_idx
        error("Unable to find block info")
    elseif length(rng) != 1
        error("Duplicated find block info")
    end
    idx = first(rng)
    return idx, block_infos[idx]
end

function diff_phase(block_infos, fit, data_idx, block_idx1, block_idx2, wrap_offset=0.5)
    idx1, block_info1 = find_block_info(block_infos, data_idx, block_idx1)
    idx2, block_info2 = find_block_info(block_infos, data_idx, block_idx2)

    # The phase (in unit of 2pi) for the two blocks are
    # `f * (t - toffset_1) + ϕf_1` and `f * (t - toffset_2) + ϕf_2`
    # Phase difference is `f * (t - toffset_2) + ϕf_2 - f * (t - toffset_1) - ϕf_1`
    # or `f * (toffset_1 - toffset_2) + ϕf_2 - ϕf_1`
    Δϕf = mod(fit.param[3] * (block_info1.time_offset - block_info2.time_offset), 1)
    Δϕf = mod(Δϕf + fit.param[3 + idx2] - fit.param[3 + idx1] + wrap_offset, 1) - wrap_offset
    vec = zeros(length(fit.param))
    vec[3] = block_info1.time_offset - block_info2.time_offset
    vec[3 + idx2] += 1
    vec[3 + idx1] += -1

    s_Δϕf = vec' * fit.covar * vec
    return Unc(Δϕf, sqrt(s_Δϕf))
end

# const block_infos_sC1, block_counts_sC1, tblocks_sC1, vblocks_sC1 =
#     collect_blocks(data_set_C1, 0.1, 0.9)
# const block_infos_sC2, block_counts_sC2, tblocks_sC2, vblocks_sC2 =
#     collect_blocks(data_set_C2, 0.1, 0.9)
# const fit_sC1 = fit_multi_phases(tblocks_sC1, vblocks_sC1)
# const fit_sC2 = fit_multi_phases(tblocks_sC2, vblocks_sC2)
# # @show fit_sC1.uncs .- fit_C1.uncs
# # @show fit_sC2.uncs .- fit_C2.uncs

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
# compare_fits(length(names_C1), block_counts_C1, block_infos_C1, fit_C1,
#              block_infos_sC1, fit_sC1)
# compare_fits(length(names_C2), block_counts_C2, block_infos_C2, fit_C2,
#              block_infos_sC2, fit_sC2)

const prefix = joinpath(@__DIR__, "data/")
const img_prefix = joinpath(@__DIR__, "imgs/")

function save_tables(prefix, names, block_infos, block_counts, fit)
    for (data_idx, name) in enumerate(names)
        name = replace(name, r"\.[^.]*$"=>"")
        open("$(prefix)_$(name).csv", "w") do io
            nblocks = block_counts[data_idx]
            print(io, "phase diff (cycle)")
            for blk_id1 in 2:nblocks
                print(io, ",$blk_id1")
            end
            println(io)
            for blk_id2 in 1:(nblocks - 1)
                print(io, blk_id2)
                for blk_id1 in 2:nblocks
                    if blk_id1 <= blk_id2
                        print(io, ",-")
                        continue
                    end
                    # Compute ϕ1 - ϕ2
                    d = diff_phase(block_infos, fit, data_idx, blk_id2, blk_id1)
                    print(io, ",$d")
                end
                println(io)
            end
        end
    end
end

save_tables("$(prefix)phase_diff", names_C1, block_infos_C1, block_counts_C1, fit_C1)
save_tables("$(prefix)phase_diff", names_C2, block_infos_C2, block_counts_C2, fit_C2)

const pairs_of_interest = [(1, 1, 2, 0),
                           # (1, 3, 4, 500),
                           (2, 1, 2, 1000),
                           (3, 1, 2, 2000),
                           (4, 1, 2, 3000),
                           (5, 1, 2, 4000),
                           (6, 1, 2, 5000),
                           (7, 1, 2, 6000),
                           (8, 1, 2, 7000),
                           (9, 1, 2, 8000),
                           (10, 1, 2, 9000)]

function model_linphase(x, p)
    return mod.(p[1] .+ p[2] .* x .+ 0.5, 1) .- 0.5
end

function save_pair_table(fname, pairs_info, block_infos_C1, block_counts_C1, fit_C1,
                         block_infos_C2, block_counts_C2, fit_C2)
    uss = Float64[]
    d_C1s = Float64[]
    d_C1_uncs = Float64[]
    d_C2s = Float64[]
    d_C2_uncs = Float64[]
    diffs = Float64[]
    diff_uncs = Float64[]
    for (data_idx, blk_id1, blk_id2, us) in pairs_info
        d_C1 = diff_phase(block_infos_C1, fit_C1, data_idx, blk_id1, blk_id2)
        d_C2 = diff_phase(block_infos_C2, fit_C2, data_idx, blk_id1, blk_id2)
        diff = d_C1 - d_C2
        diff = Unc(mod(diff.a + 0.5, 1) - 0.5, diff.s)
        push!(uss, us)
        push!(d_C1s, d_C1.a)
        push!(d_C1_uncs, d_C1.s)
        push!(d_C2s, d_C2.a)
        push!(d_C2_uncs, d_C2.s)
        push!(diffs, diff.a)
        push!(diff_uncs, diff.s)
    end

    fit_C1 = fit_data(model_linphase, uss, d_C1s, d_C1_uncs, [0.0, 0.00005], plotx=false)
    fit_C2 = fit_data(model_linphase, uss, d_C2s, d_C2_uncs, [0.0, 0.00005], plotx=false)

    for i in 1:length(uss)
        d_C1s[i] = mod(d_C1s[i] - fit_C1.param[2] * uss[i] + 0.5, 1) - 0.5
        d_C2s[i] = mod(d_C2s[i] - fit_C2.param[2] * uss[i] + 0.5, 1) - 0.5
    end

    open(fname, "w") do io
        println(io, "time (us),phase ind1 (cycle),phase ind2 (cycle),phase diff (cycle)")
        for (us, d_C1, d_C1_unc, d_C2, d_C2_unc, diff, diff_unc) in zip(uss, d_C1s, d_C1_uncs, d_C2s, d_C2_uncs, diffs, diff_uncs)
            println(io, "$us,$(Unc(d_C1, d_C1_unc)),$(Unc(d_C2, d_C2_unc)),$(Unc(diff, diff_unc))")
        end
    end
    return uss, d_C1s, d_C1_uncs, d_C2s, d_C2_uncs, diffs, diff_uncs
end

const uss, d_C1s, d_C1_uncs, d_C2s, d_C2_uncs, diffs, diff_uncs =
    save_pair_table("$(prefix)pulse_phase_20221108.csv",
                    pairs_of_interest, block_infos_C1, block_counts_C1, fit_C1,
                    block_infos_C2, block_counts_C2, fit_C2)


figure()
errorbar(uss, d_C1s, d_C1_uncs, label="ind1")
errorbar(uss, d_C2s, d_C2_uncs, label="ind2")
errorbar(uss, diffs, diff_uncs, label="diff")
grid()
legend(ncol=3, fontsize=12)
xlabel("Wait time (\$\\mu s\$)")
ylabel("Phase (cycle)")
tight_layout()
NaCsPlot.maybe_save("$(img_prefix)pulse_phase_20221108")

NaCsPlot.maybe_show()

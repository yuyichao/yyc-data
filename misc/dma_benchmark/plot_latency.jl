#!/usr/bin/julia

include("utils.jl")

push!(Base.LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using YAML
using NaCsPlot
using PyPlot

const all_lines = Dict{Any,LineGroup{HistCollector}}()

for f in readdir(joinpath(@__DIR__, "data/latency"), join=true)
    for d in YAML.load_file(f)
        add_nbuff_item!(all_lines, d)
    end
end

const prefix = joinpath(@__DIR__, "imgs/dmatest")

const cm = get(matplotlib.colormaps, "viridis")

figure(figsize=(6.4 * 4, 4.8 * 4))
for (i, k) in enumerate(("HP DDR", "HP DDR_WC", "HP OCM", "HP OCM_WC",
                         "ACP DDR", "ACP DDR_WC", "ACP OCM", "ACP OCM_WC",
                         "ACP_L2 DDR", "ACP_L2 DDR_WC", "ACP_L2 OCM", "ACP_L2 OCM_WC",
                         "ACP_L1 DDR", "ACP_L1 DDR_WC", "ACP_L1 OCM", "ACP_L1 OCM_WC",))
    subplot(4, 4, i)
    cs = all_lines["dma_latency"].lines[k].collectors
    nbuffs = sort!(collect(keys(cs)))
    max_nbuff = nbuffs[end]
    n_nbuff = length(nbuffs)
    min_val = typemax(Int)
    max_val = typemin(Int)
    for (i, nbuff) in enumerate(nbuffs)
        freq = cs[nbuff].freq
        vals = sort!(collect(keys(freq)))
        min_val = min(vals[1], min_val)
        max_val = max(vals[end], max_val)
        freqs = [freq[val] for val in vals]
        bar(vals .+ ((i - 1) / (n_nbuff - 1) - 0.5) * 0.9, freqs ./ sum(freqs), width=1.8 / n_nbuff, color=cm(nbuff / max_nbuff), alpha=0.8)
    end
    grid()
    xlim([min_val - 0.49, max_val + 0.49])
    ylim([0, ylim()[2]])
    title(endswith(k, "_WC") ? (k[1:end - 3] * " CLR") : k)
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_dma_latency_hist")

figure(figsize=(6.4 * 2, 4.8))
subplot(1, 2, 1)
for k in ("HP DDR", "HP DDR_WC", "ACP DDR", "ACP DDR_WC",
          "ACP_L2 DDR", "ACP_L2 DDR_WC", "ACP_L1 DDR", "ACP_L1 DDR_WC")
    line = all_lines["dma_latency"].lines[k]
    s, v, u = get_data(line)
    errorbar(s, v, u, label=endswith(k, "_WC") ? (k[1:end - 3] * " CLR") : k)
end
grid()
xlim([0, 70])
xlabel("Buffer Number")
ylabel("Latency (cycle)")
legend(fontsize=11, ncol=2)
title("DMA Latency (DDR)")

subplot(1, 2, 2)
for k in ("HP OCM", "HP OCM_WC", "ACP OCM", "ACP OCM_WC",
          "ACP_L2 OCM", "ACP_L2 OCM_WC", "ACP_L1 OCM", "ACP_L1 OCM_WC")
    line = all_lines["dma_latency"].lines[k]
    s, v, u = get_data(line)
    errorbar(s, v, u, label=endswith(k, "_WC") ? (k[1:end - 3] * " CLR") : k)
end
grid()
xlim([0, 70])
xlabel("Buffer Number")
ylabel("Latency (cycle)")
legend(fontsize=11, ncol=2)
title("DMA Latency (OCM)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_dma_latency")

NaCsPlot.maybe_show()

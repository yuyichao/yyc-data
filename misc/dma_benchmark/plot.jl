#!/usr/bin/julia

include("utils.jl")

push!(Base.LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using YAML
using NaCsPlot
using PyPlot

const all_lines = Dict{Any,LineGroup{Collector}}()

for f in readdir(joinpath(@__DIR__, "data/throughput"), join=true)
    for d in YAML.load_file(f)
        add_item!(all_lines, d)
    end
end

const prefix = joinpath(@__DIR__, "imgs/dmatest")

figure()
for k in ("DDR", "DDR_WC", "OCM", "OCM_WC")
    line = all_lines["memcpy"].lines[k]
    s, v, u = get_data(line)
    errorbar(s, s ./ v, u ./ v.^2 .* s, label=k)
end
grid()
xscale("log")
xlim([128, 65536])
xticks([256, 1024, 4096, 16384, 65536],
       ["256", "1k", "4k", "16k", "64k"])
xticks([512, 2048, 8192, 32768], minor=true)
xlabel("Block size (B)")
ylabel("Throughput (B/ns)")
legend(fontsize=11, ncol=2)
title("memcpy")
NaCsPlot.maybe_save("$(prefix)_memcpy")

figure()
for k in ("DDR", "DDR_WC", "OCM", "OCM_WC")
    line = all_lines["rand_fill"].lines[k]
    s, v, u = get_data(line)
    errorbar(s, s ./ v, u ./ v.^2 .* s, label=k)
end
grid()
xscale("log")
xlim([128, 65536])
xticks([256, 1024, 4096, 16384, 65536],
       ["256", "1k", "4k", "16k", "64k"])
xticks([512, 2048, 8192, 32768], minor=true)
xlabel("Block size (B)")
ylabel("Throughput (B/ns)")
legend(fontsize=11, ncol=2)
title("rand fill")
NaCsPlot.maybe_save("$(prefix)_rand")

figure(figsize=(6.4 * 2, 4.8))
subplot(1, 2, 1)
for k in ("HP DDR_WC", "HP OCM_WC",
          "ACP DDR_WC", "ACP OCM_WC",
          "ACP_L2 DDR_WC", "ACP_L2 OCM_WC",
          "ACP_L1 DDR", "ACP_L1 DDR_WC", "ACP_L1 OCM", "ACP_L1 OCM_WC")
    line = all_lines["flush"].lines[k]
    s, v, u = get_data(line)
    if k == "ACP_L2 DDR_WC" || k == "ACP_L1 DDR_WC"
        errorbar(s, v, u, label=k, ls="--")
    else
        errorbar(s, v, u, label=k)
    end
end
grid()
xscale("log")
xlim([128, 65536])
xticks([256, 1024, 4096, 16384, 65536],
       ["256", "1k", "4k", "16k", "64k"])
xticks([512, 2048, 8192, 32768], minor=true)
xlabel("Block size (B)")
ylim([8, 10])
ylabel("Latency (ns)")
legend(fontsize=11, ncol=2)
title("DMA prep (fast)")

subplot(1, 2, 2)
for k in ("HP DDR", "HP OCM",
          "ACP DDR", "ACP OCM",
          "ACP_L2 DDR", "ACP_L2 OCM")
    line = all_lines["flush"].lines[k]
    s, v, u = get_data(line)
    errorbar(s, s ./ v, u ./ v.^2 .* s, label=k)
end
grid()
xscale("log")
xlim([128, 65536])
xticks([256, 1024, 4096, 16384, 65536],
       ["256", "1k", "4k", "16k", "64k"])
xticks([512, 2048, 8192, 32768], minor=true)
xlabel("Block size (B)")
ylabel("Throughput (B/ns)")
legend(fontsize=11, ncol=2)
title("DMA prep (slow)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_flush")

figure(figsize=(6.4 * 2, 4.8 * 4))
for i in 1:4
    subplot(4, 2, i * 2 - 1)
    for k in ("HP DDR_WC", "HP OCM_WC",
              "ACP DDR_WC", "ACP OCM_WC",
              "ACP_L2 DDR_WC", "ACP_L2 OCM_WC",
              "ACP_L1 DDR", "ACP_L1 DDR_WC", "ACP_L1 OCM", "ACP_L1 OCM_WC")
        line = all_lines[("dma_pipe", i)].lines[k]
        s, v, u = get_data(line)
        if k == "ACP_L2 DDR_WC" || k == "ACP_L1 DDR_WC"
            errorbar(s, s ./ v, u ./ v.^2 .* s, label=k, ls="--")
        else
            errorbar(s, s ./ v, u ./ v.^2 .* s, label=k)
        end
    end
    axhline(1.6, color="red", ls="-.", alpha=0.6)
    grid()
    xscale("log")
    xlim([128, 65536])
    ylim([0, 1.63])
    xticks([256, 1024, 4096, 16384, 65536],
           ["256", "1k", "4k", "16k", "64k"])
    xticks([512, 2048, 8192, 32768], minor=true)
    xlabel("Block size (B)")
    ylabel("Throughput (B/ns)")
    legend(fontsize=11, ncol=2, loc="lower right")
    title("Pipe $i blocks (fast prep)")

    subplot(4, 2, i * 2)
    for k in ("HP DDR", "HP OCM",
              "ACP DDR", "ACP OCM",
              "ACP_L2 DDR", "ACP_L2 OCM")
        line = all_lines[("dma_pipe", i)].lines[k]
        s, v, u = get_data(line)
        errorbar(s, s ./ v, u ./ v.^2 .* s, label=k)
    end
    grid()
    xscale("log")
    xlim([128, 65536])
    ylim([0, 0.55])
    xticks([256, 1024, 4096, 16384, 65536],
           ["256", "1k", "4k", "16k", "64k"])
    xticks([512, 2048, 8192, 32768], minor=true)
    xlabel("Block size (B)")
    ylabel("Throughput (B/ns)")
    legend(fontsize=11, ncol=2, loc="lower right")
    title("Pipe $i blocks (slow prep)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_dma_pipe")

figure(figsize=(6.4 * 2, 4.8 * 4))
for i in 1:4
    subplot(4, 2, i * 2 - 1)
    for k in ("HP DDR", "HP DDR_WC", "HP OCM", "HP OCM_WC",
              "ACP DDR", "ACP DDR_WC", "ACP OCM", "ACP OCM_WC",
              "ACP_L2 OCM", "ACP_L2 OCM_WC")
        line = all_lines[("dma_only", i)].lines[k]
        s, v, u = get_data(line)
        if k == "ACP_L2 DDR_WC" || k == "ACP_L1 DDR_WC"
            errorbar(s, s ./ v, u ./ v.^2 .* s, label=k, ls="--")
        else
            errorbar(s, s ./ v, u ./ v.^2 .* s, label=k)
        end
    end
    axhline(1.6, color="red", ls="-.", alpha=0.6)
    grid()
    xscale("log")
    xlim([128, 65536])
    ylim([0, 1.63])
    xticks([256, 1024, 4096, 16384, 65536],
           ["256", "1k", "4k", "16k", "64k"])
    xticks([512, 2048, 8192, 32768], minor=true)
    xlabel("Block size (B)")
    ylabel("Throughput (B/ns)")
    legend(fontsize=11, ncol=2, loc="lower right")
    title("DMA $i blocks (full speed)")

    subplot(4, 2, i * 2)
    for k in ("ACP_L2 DDR", "ACP_L2 DDR_WC",
              "ACP_L1 DDR", "ACP_L1 DDR_WC", "ACP_L1 OCM", "ACP_L1 OCM_WC")
        line = all_lines[("dma_only", i)].lines[k]
        s, v, u = get_data(line)
        if k == "ACP_L2 DDR_WC" || k == "ACP_L1 DDR_WC"
            errorbar(s, s ./ v, u ./ v.^2 .* s, label=k, ls="--")
        else
            errorbar(s, s ./ v, u ./ v.^2 .* s, label=k)
        end
    end
    axhline(1.6, color="red", ls="-.", alpha=0.6)
    grid()
    xscale("log")
    xlim([128, 65536])
    ylim([0, 1.63])
    xticks([256, 1024, 4096, 16384, 65536],
           ["256", "1k", "4k", "16k", "64k"])
    xticks([512, 2048, 8192, 32768], minor=true)
    xlabel("Block size (B)")
    ylabel("Throughput (B/ns)")
    legend(fontsize=11, ncol=2, loc="lower right")
    title("DMA $i blocks (slower)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_dma")

NaCsPlot.maybe_show()

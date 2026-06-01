#!/usr/bin/julia

include("utils.jl")

push!(Base.LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using YAML
using NaCsPlot
using PyPlot

for f in readdir(joinpath(@__DIR__, "data"), join=true)
    for d in YAML.load_file(f)
        add_item!(d)
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
xlabel("Block size (B)")
ylabel("Throughput (B/ns)")
legend(fontsize=12, ncol=2)
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
xlabel("Block size (B)")
ylabel("Throughput (B/ns)")
legend(fontsize=12, ncol=2)
NaCsPlot.maybe_save("$(prefix)_rand")

NaCsPlot.maybe_show()


# ("dma_pipe", 1)
# ("dma_pipe", 2)
# ("dma_pipe", 3)
# ("dma_pipe", 4)

# ("dma_only", 1)
# ("dma_only", 2)
# ("dma_only", 3)
# ("dma_only", 4)

# "flush"

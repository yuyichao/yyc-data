#!/usr/bin/julia

include("process_data.jl")

using NaCsPlot
using PyPlot

const prefix = joinpath(@__DIR__, "imgs/freqs_20230915_2")

const fits = [@time(fit_freqs(read_bin_compressed(joinpath(@__DIR__, "data/$(i).zst")),
                              6000)) for i in 60:78]

figure(figsize=[6.4 * 4, 4.8 * 3])
for i in 2:19
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([0, 5])
grid()
xlabel("Time (ms)")
ylabel("Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_full1")

figure(figsize=[6.4 * 4, 4.8])
for i in 1:1
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([0, 5])
grid()
xlabel("Time (ms)")
ylabel("Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_full2")

NaCsPlot.maybe_show()

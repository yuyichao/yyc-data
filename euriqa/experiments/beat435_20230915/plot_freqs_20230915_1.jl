#!/usr/bin/julia

include("process_data.jl")

using NaCsPlot
using PyPlot

const prefix = joinpath(@__DIR__, "imgs/freqs_20230915_1")

const fits = [@time(fit_freqs(read_bin_compressed(joinpath(@__DIR__, "data/$(i).zst")),
                              6000)) for i in 79:91]

figure(figsize=[6.4 * 4, 4.8 * 3])
subplot(3, 1, 1)
for i in 1:3:13
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([0, 5])
grid()
# xlabel("Time (ms)")
ylabel("Frequency (MHz)")

subplot(3, 1, 2)
for i in 2:3:13
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([0, 5])
grid()
# xlabel("Time (ms)")
ylabel("Frequency (MHz)")

subplot(3, 1, 3)
for i in 3:3:13
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([0, 5])
grid()
xlabel("Time (ms)")
ylabel("Frequency (MHz)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_full")

figure(figsize=[6.4 * 4, 4.8 * 3])
subplot(3, 1, 1)
for i in 1:3:13
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([1.5, 3.5])
grid()
# xlabel("Time (ms)")
ylabel("Frequency (MHz)")

subplot(3, 1, 2)
for i in 2:3:13
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([1.5, 3.5])
grid()
# xlabel("Time (ms)")
ylabel("Frequency (MHz)")

subplot(3, 1, 3)
for i in 3:3:13
    errorbar(fits[i].ts .* 1e3, fits[i].freqs ./ 1e6, fits[i].freqs_unc ./ 1e6)
end
xlim([1.5, 3.5])
grid()
xlabel("Time (ms)")
ylabel("Frequency (MHz)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_zoomin")

NaCsPlot.maybe_show()

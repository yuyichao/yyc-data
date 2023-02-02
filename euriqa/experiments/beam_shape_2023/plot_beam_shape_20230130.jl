#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data
using NaCsPlot
using PyPlot
import NaCsCalc.Format: Unc, Sci

function load_count(pos, vertical)
    name = joinpath(@__DIR__, "data",
                    "counts_$(pos)_$(vertical ? 'v' : 'h')_20230130.csv")
    return readdlm(name, ',', Float64, skipstart=1)
end

const zpositions = [10.0, 10.5, 11.0, 11.5, 12.075, 12.5]

const hdata = [load_count(pos, false) for pos in zpositions]
const vdata = [load_count(pos, true) for pos in zpositions]

const prefix = joinpath(@__DIR__, "imgs", "beam_shape_20230130")

function model(x, p)
    return p[1] .+ p[2] .* exp.(.-2 .* (x .- p[3]).^2 ./ p[4]^2)
end

function fit_count(data)
    return fit_data(model, data[:, 1], data[:, 2],
                    [0, 20, (data[1, 1] + data[end, 1]) / 2,
                     (data[end, 1] - data[1, 1]) / 2])
end

const hfit = [fit_count(data) for data in hdata]
const vfit = [fit_count(data) for data in vdata]

figure(figsize=[6.4 * 2, 4.8])
subplot(1, 2, 1)
for i in 1:length(zpositions)
    data = hdata[i]
    fit = hfit[i]
    plot(data[:, 1], data[:, 2], "C$(i - 1)o")
    plot(fit.plotx, fit.ploty, "C$(i - 1)")
end
xlim([6.83, 7.23])
ylim([1.25, 54])
text(6.83, 50, "\$w_{$(zpositions[1])}=$(hfit[1].uncs[4] * 1000)\\ \\mu m\$",
     color="C0", size=16)
text(6.83, 46, "\$w_{$(zpositions[2])}=$(hfit[2].uncs[4] * 1000)\\ \\mu m\$",
     color="C1", size=16)
text(6.83, 42, "\$w_{$(zpositions[3])}=$(hfit[3].uncs[4] * 1000)\\ \\mu m\$",
     color="C2", size=16)
text(6.83, 38, "\$w_{$(zpositions[4])}=$(hfit[4].uncs[4] * 1000)\\ \\mu m\$",
     color="C3", size=16)
text(6.83, 34, "\$w_{$(zpositions[5])}=$(hfit[5].uncs[4] * 1000)\\ \\mu m\$",
     color="C4", size=16)
text(6.83, 30, "\$w_{$(zpositions[6])}=$(hfit[6].uncs[4] * 1000)\\ \\mu m\$",
     color="C5", size=16)
grid()
title("Horizontal Profile")
xlabel("Position (mm)")
ylabel("Counts")

subplot(1, 2, 2)
for i in 1:length(zpositions)
    data = vdata[i]
    fit = vfit[i]
    plot(data[:, 1], data[:, 2], "C$(i - 1)o")
    plot(fit.plotx, fit.ploty, "C$(i - 1)")
end
xlim([6.425, 6.4885])
ylim([3, 54])
text(6.425, 50, "\$w_{$(zpositions[1])}=$(vfit[1].uncs[4] * 1000)\\ \\mu m\$",
     color="C0", size=16)
text(6.425, 46, "\$w_{$(zpositions[2])}=$(vfit[2].uncs[4] * 1000)\\ \\mu m\$",
     color="C1", size=16)
text(6.425, 42, "\$w_{$(zpositions[3])}=$(vfit[3].uncs[4] * 1000)\\ \\mu m\$",
     color="C2", size=16)
text(6.425, 38, "\$w_{$(zpositions[4])}=$(vfit[4].uncs[4] * 1000)\\ \\mu m\$",
     color="C3", size=16)
text(6.425, 34, "\$w_{$(zpositions[5])}=$(vfit[5].uncs[4] * 1000)\\ \\mu m\$",
     color="C4", size=16)
text(6.425, 30, "\$w_{$(zpositions[6])}=$(vfit[6].uncs[4] * 1000)\\ \\mu m\$",
     color="C5", size=16)
grid()
title("Vertical Profile")
xlabel("Position (mm)")
ylabel("Counts")

tight_layout()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

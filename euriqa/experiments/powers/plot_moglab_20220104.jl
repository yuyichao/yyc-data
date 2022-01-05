#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
import NaCsCalc.Format: Unc, Sci
using NaCsData.Fitting: fit_data
using Printf

const iname369 = joinpath(@__DIR__, "data", "moglab_369_20220104.csv")
const data369 = readdlm(iname369, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "moglab_20220104")

model_lin(x, p) = p[1] .+ p[2] .* x

const fit369_idxs = data369[:, 2] .> 2
const fit369 = fit_data(model_lin, data369[fit369_idxs, 1], data369[fit369_idxs, 2],
                        [-60, 1.0], plot_lo=60)

figure()
plot(data369[:, 1], data369[:, 2], "C0o-")
plot(fit369.plotx, fit369.ploty, "C2--")
text(62.8, 4, "$(@sprintf("%.2f", fit369.param[2])) mW/mA", color="C2",
     rotation=atand(fit369.param[2]), transform_rotates_text=true)
ylim([0, 12.9])
grid()
xlabel("Current (mA)")
ylabel("Output Power (mW)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_369")

NaCsPlot.maybe_show()

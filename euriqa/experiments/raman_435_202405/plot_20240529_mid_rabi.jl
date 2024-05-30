#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_loading
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using Printf
using SpecialFunctions

const name_rabi, (_, _, data_rabi,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056939-Raman435MidCircuitScanPreRabi.h5"),
                            thresh=2)

const prefix = joinpath(@__DIR__, "imgs", "data_20240529_mid_rabi")

function model(x, p)
    return p[1] .- p[2] / 2 .* cos.(2π .* x)
end

const fit_rabi = fit_loading(model, data_rabi, [0.5, 1])
@show fit_rabi.uncs

figure()
NaCsPlot.plot_loading_data(data_rabi, xscale=360, fmt="C0o")
plot(fit_rabi.plotx .* 360, fit_rabi.ploty, color="C0")
grid()
ylim([0, 1])
xlim([0, 360])
xlabel("Single qubit rotation angle (\$^\\circ\$)")
ylabel("\$P_1\$")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

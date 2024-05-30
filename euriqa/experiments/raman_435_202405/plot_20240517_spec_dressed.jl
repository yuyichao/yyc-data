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

const name_right, (_, data_right,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000054729-Raman435ScanSpec.h5"),
                            thresh=2)
const name_left, (_, data_left,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000054730-Raman435ScanSpec.h5"),
                            thresh=2)
const name_center, (_, data_center,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000054731-Raman435ScanSpec.h5"),
                            thresh=2)

const prefix = joinpath(@__DIR__, "imgs", "data_20240524_spec_dressed")

function rabiLine(det, t, Omega)
    Omega2 = Omega^2
    OmegaG2 = det^2 + Omega2
    return Omega2 / OmegaG2 * sin(√(OmegaG2) * t / 2)^2
end

function gen_model(t)
    return (x, p) -> p[1] .+ p[2] .* rabiLine.(2π .* (x .- p[3]), t, p[4])
end

const fit_right = fit_loading(gen_model(10e-6), data_right,
                              [1, -0.9, 200e3, π / 10e-6])
const fit_left = fit_loading(gen_model(10e-6), data_left,
                              [1, -0.9, -200e3, π / 10e-6])
const fit_center = fit_loading(gen_model(10e-6), data_center,
                              [1, -0.9, 0, π / 10e-6])
@show fit_right.uncs
@show fit_left.uncs
@show fit_center.uncs

const fcenter = fit_center.param[3]

const data_left_shift = NaCsData.map_params((i, v)->v - fcenter, data_left)
const data_right_shift = NaCsData.map_params((i, v)->v - fcenter, data_right)
const data_center_shift = NaCsData.map_params((i, v)->v - fcenter, data_center)

figure()
NaCsPlot.plot_loading_data(data_left_shift, -1, yoffset=1, xscale=1e-3, fmt="C0o",
                           label="Dress \$+y\$")
plot((fit_left.plotx .- fcenter) ./ 1000, 1 .- fit_left.ploty, color="C0")
NaCsPlot.plot_loading_data(data_center_shift, -1, yoffset=1, xscale=1e-3, fmt="C1o",
                           label="No dressing")
plot((fit_center.plotx .- fcenter) ./ 1000, 1 .- fit_center.ploty, color="C1")
NaCsPlot.plot_loading_data(data_right_shift, -1, yoffset=1, xscale=1e-3, fmt="C2o",
                           label="Dress \$-y\$")
plot((fit_right.plotx .- fcenter) ./ 1000, 1 .- fit_right.ploty, color="C2")
grid()
ylim([0, 1.06])
xlim([-350, 350])
legend(ncol=3, fontsize=11, loc="upper center")
xlabel("Detuning from resonance (kHz)")
ylabel("\$P_{D_{3/2}}\$")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

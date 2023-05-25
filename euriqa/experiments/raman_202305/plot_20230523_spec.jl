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

const name_spec, (data_spec,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000018928-RamanSpecGateScan.h5"),
                            thresh=31)

const prefix = joinpath(@__DIR__, "imgs", "data_20230523_spec")

function rabiLine(det, t, Omega)
    Omega2 = Omega^2
    OmegaG2 = det^2 + Omega2
    return Omega2 / OmegaG2 * sin(√(OmegaG2) * t / 2)^2
end

function gen_model(t)
    return (x, p) -> p[1] .+ p[2] .* rabiLine.(2π .* (x .- p[3]), t, p[4])
end

const fit_spec = fit_loading(gen_model(90e-6), data_spec,
                             [0.02, 0.9, 209.18e6, π / 90e-6])
@show fit_spec.uncs

const data_spec_shift = NaCsData.map_params((i, v)->v - fit_spec.param[3], data_spec)


figure()
NaCsPlot.plot_loading_data(data_spec_shift, xscale=1e-3, fmt="C0o")
plot((fit_spec.plotx .- fit_spec.param[3]) / 1000, fit_spec.ploty, color="C0")
grid()
ylim([0, 0.85])
xlabel("Detuning from resonance (kHz)")
ylabel("\$P_{1}\$")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

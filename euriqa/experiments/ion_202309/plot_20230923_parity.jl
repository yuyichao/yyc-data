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

const param_name, (params, logicals) =
    NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data/000030222-DriveTwoQubitScan.h5"))

function even_selector(logicals)
    @assert size(logicals) == (1, 2)
    return [1, logicals[1, 1] == logicals[1, 2]]
end

const data = NaCsData.select_count(params, logicals, even_selector)

const prefix = joinpath(@__DIR__, "imgs", "data_20230923_parity_30222")

function sin2(det, t, Omega)
    Omega2 = Omega^2
    OmegaG2 = det^2 + Omega2
    return Omega2 / OmegaG2 * sin(√(OmegaG2) * t / 2)^2
end

function model(x, p)
    return p[1] .+ p[2] .* sin.(4π .* (x .- p[3]))
end

const fit = fit_loading(model, data, [0.55, 0.4, 0.3])
@show fit.uncs

figure(figsize=[6.4 * 2, 4.8])
NaCsPlot.plot_loading_data(data, 2, yoffset=-1, xscale=2π, fmt="C0o")
plot(fit.plotx .* 2π, fit.ploty .* 2 .- 1, color="C0")
grid()
xlim([0, 2π])
# ylim([-1, 1])
xlabel("\$\\phi\$ (rad)")
ylabel("\$P_{11,00}-P_{10,01}\$")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

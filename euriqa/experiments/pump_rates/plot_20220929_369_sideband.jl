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

const name_369_1_1, data_369_1_1 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000006041-F1OPRateScan.h5"))
const name_369_2, data_369_2 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000006042-F1OPRateScan.h5"))
const name_369_3_1, data_369_3_1 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000006043-F1OPRateScan.h5"))
const name_369_3_2, data_369_3_2 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000006044-F1OPRateScan.h5"))
const name_369_4, data_369_4 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000006045-F1OPRateScan.h5"))
const name_369_1_2, data_369_1_2 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000006046-F1OPRateScan.h5"))

const data_369_1 = [data_369_1_1; data_369_1_2]
const data_369_3 = [data_369_3_1; data_369_3_2]

const prefix = joinpath(@__DIR__, "imgs", "data_20220929_369_sideband")

model_exp(xs, p) = p[1] .+ p[2] .* exp.(.-xs ./ p[3])

const fit_369_1 = fit_loading(model_exp, data_369_1, [0.02, 2, 0.6e-6])
const fit_369_2 = fit_loading(model_exp, data_369_2, [0.02, 2, 0.7e-6])
const fit_369_3 = fit_loading(model_exp, data_369_3, [0.02, 2, 0.8e-6])
const fit_369_4 = fit_loading(model_exp, data_369_4, [0.02, 2, 0.9e-6])

const powers_dBm = [-0.26, -2.55, -3.49, -4.43] .+ 30.23
const powers_W = 10 .^ (powers_dBm ./ 10) ./ 1000
const powers_name = [@sprintf("%.2f W", p) for p in powers_W]

figure()
NaCsPlot.plot_loading_data(data_369_1, xscale=1e6, fmt="C0o", label=powers_name[1])
plot(fit_369_1.plotx .* 1e6, fit_369_1.ploty, color="C0")
NaCsPlot.plot_loading_data(data_369_2, xscale=1e6, fmt="C1o", label=powers_name[2])
plot(fit_369_2.plotx .* 1e6, fit_369_2.ploty, color="C1")
NaCsPlot.plot_loading_data(data_369_3, xscale=1e6, fmt="C2o", label=powers_name[3])
plot(fit_369_3.plotx .* 1e6, fit_369_3.ploty, color="C2")
NaCsPlot.plot_loading_data(data_369_4, xscale=1e6, fmt="C3o", label=powers_name[4])
plot(fit_369_4.plotx .* 1e6, fit_369_4.ploty, color="C3")
legend(ncol=2, fontsize=12)
text(2.3, 0.25, "\$\\tau_{\\mathrm{$(powers_name[1])}}=$(fit_369_1.uncs[3] * 1e6) \\mu s\$", color="C0")
text(2.3, 0.34, "\$\\tau_{\\mathrm{$(powers_name[2])}}=$(fit_369_2.uncs[3] * 1e6) \\mu s\$", color="C1")
text(2.3, 0.43, "\$\\tau_{\\mathrm{$(powers_name[3])}}=$(fit_369_3.uncs[3] * 1e6) \\mu s\$", color="C2")
text(2.3, 0.54, "\$\\tau_{\\mathrm{$(powers_name[4])}}=$(fit_369_4.uncs[3] * 1e6) \\mu s\$", color="C3")
grid()
xlim([0.5, 5])
ylim([0, 1])
xlabel("Pump Time (\$\\mu s\$)")
ylabel("\$P_{detect}\$")
NaCsPlot.maybe_save("$(prefix)_rate")

const taus_uncs = [fit_369_1.uncs[3], fit_369_2.uncs[3],
                   fit_369_3.uncs[3], fit_369_4.uncs[3]] .* 1e6
const rates_uncs = 1 ./ taus_uncs
const rates = [r.a for r in rates_uncs]
const rates_unc = [r.s for r in rates_uncs]

model_eom(xs, p) = p[1] .* besselj1.(sqrt.(xs ./ p[2])).^2

const fit_eom = fit_data(model_eom, powers_W, rates, rates_unc, [7.0, 0.7])

figure()
errorbar(powers_W, rates, rates_unc, fmt="C0o")
plot(fit_eom.plotx, fit_eom.ploty, color="C0")
text(0.63, 1.05, "\$P_{\\mathrm{1rad}}=$(fit_eom.uncs[2]) W\$", color="C0")
xlabel("RF Power (W)")
ylabel("Pump Rate (\$\\mu s^{-1}\$)")
grid()
NaCsPlot.maybe_save("$(prefix)_power")


NaCsPlot.maybe_show()

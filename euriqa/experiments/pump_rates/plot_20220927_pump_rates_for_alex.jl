#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_loading
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const name_369, data_369 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000003835-F1OPRateScan.h5"))
const name_935_1, data_935_1 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000003174-DF2OPRateScan.h5"))
const name_935_2, data_935_2 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000003176-DF2OPRateScan.h5"))
const name_935_3, data_935_3 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000003177-DF2OPRateScan.h5"))
const name_935fix, data_935fix =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000005317-DF2OPRateScan.h5"))
const name_depump0, data_depump0 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000003179-MFDepumpRateScan.h5"))
const name_depump1, data_depump1 =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000003668-MFDepumpRateScan.h5"))

const data_935 = [data_935_1; data_935_2; data_935_3;]

const prefix = joinpath(@__DIR__, "imgs", "data_20220927_pump_rates_for_alex")

model_exp(xs, p) = p[1] .+ p[2] .* exp.(.-xs ./ p[3])

function gen_exp1_mod(p0)
    (xs, p) -> p0 .+ p[1] .* exp.(.-xs ./ p[2])
end

const fit_369 = fit_loading(model_exp, data_369[3:20], [0.05, 0.95, 1e-6])
const fit_935 = fit_loading(model_exp, data_935, [0.7, -0.7, 10e-3])
const fit_935fix = fit_loading(gen_exp1_mod(0.60), data_935fix, [-0.55, 0.4e-6])
const fit_depump0 = fit_loading(model_exp, data_depump0, [0.7, -0.7, 200e-3])
const fit_depump1 = fit_loading(gen_exp1_mod(0.842), data_depump1, [-0.7, 2])

figure()
NaCsPlot.plot_loading_data(data_369, xscale=1e6, fmt="C0o")
plot(fit_369.plotx .* 1e6, fit_369.ploty, color="C0")
text(1.5, 0.47, "\$\\tau=$(fit_369.uncs[3] * 1e6) \\mu s\$", color="C0")
grid()
xlim([0, 5])
ylim([0, 1])
xlabel("Pump Time (\$\\mu s\$)")
ylabel("\$P_{detect}\$")
NaCsPlot.maybe_save("$(prefix)_369")

figure()
NaCsPlot.plot_loading_data(data_935, xscale=1e3, fmt="C0o")
plot(fit_935.plotx .* 1e3, fit_935.ploty, color="C0")
text(3, 0.2, "\$\\tau=$(fit_935.uncs[3] * 1e3) ms\$", color="C0")
grid()
xlim([0, 10])
xlabel("Pump Time (\$ms\$)")
ylabel("\$P_{detect}\$")
NaCsPlot.maybe_save("$(prefix)_935")

figure()
NaCsPlot.plot_loading_data(data_935fix, xscale=1e6, fmt="C0o")
plot(fit_935fix.plotx .* 1e6, fit_935fix.ploty, color="C0")
text(0.8, 0.22, "\$\\tau=$(fit_935fix.uncs[2] * 1e6) \\mu s\$", color="C0")
grid()
xlim([0, 2.5])
xlabel("Pump Time (\$\\mu s\$)")
ylabel("\$P_{detect}\$")
NaCsPlot.maybe_save("$(prefix)_935fix")

figure()
NaCsPlot.plot_loading_data(data_depump0, fmt="C0o")
plot(fit_depump0.plotx, fit_depump0.ploty, color="C0")
NaCsPlot.plot_loading_data(data_depump1, fmt="C1o")
plot(fit_depump1.plotx, fit_depump1.ploty, color="C1")
text(0.1, 0.4, "\$\\tau_{before}=$(fit_depump0.uncs[3] * 1e3) ms\$", color="C0")
text(0.1, 0.25, "\$\\tau_{after}=$(fit_depump1.uncs[2]) s\$", color="C1")
grid()
# xlim([0, 2.5])
xlabel("Pump Time (\$s\$)")
ylabel("\$P_{detect}\$")
NaCsPlot.maybe_save("$(prefix)_depump")

NaCsPlot.maybe_show()

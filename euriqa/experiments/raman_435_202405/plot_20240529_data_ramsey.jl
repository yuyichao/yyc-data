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

const xc = 0.04

const detect_crosstalk = [1 xc 0
                          xc 1 xc
                          0 xc 1]
const detect_disentangle = inv(detect_crosstalk)

function count_processor!(counts)
    counts .= detect_disentangle * counts
end

const name_lo130, (_, data_lo130,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056795-Raman435MidCircuitScanRamsey.h5"),
                            thresh=2, count_processor=count_processor!)
const name_lo300, (_, data_lo300,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056821-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_lo600, (_, data_lo600,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056824-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_lo1200, (_, data_lo1200,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056825-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_lo2400, (_, data_lo2400,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056830-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_lo4800, (_, data_lo4800,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056831-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_lo9600, (_, data_lo9600,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056832-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_lo14400, (_, data_lo14400,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056838-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_lo19200, (_, data_lo19200,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056836-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)

const name_hi130, (_, data_hi130, aux_hi130) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056796-Raman435MidCircuitScanRamsey.h5"),
                            thresh=2, count_processor=count_processor!)
const name_hi300, (_, data_hi300, aux_hi300) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056822-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_hi600, (_, data_hi600, aux_hi600) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056823-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_hi1200, (_, data_hi1200, aux_hi1200) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056826-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_hi2400, (_, data_hi2400, aux_hi2400) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056828-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)
const name_hi9600, (_, data_hi9600, aux_hi9600) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__,
                                     "data/000056835-Raman435MidCircuitScanRamsey.h5"),
                            thresh=8, count_processor=count_processor!)

const prefix = joinpath(@__DIR__, "imgs", "data_20240529_data_ramsey")

function model(x, p)
    return p[1] .+ p[2] / 2 .* cos.(2π .* (x .- p[3]))
end

const fit_lo130 = fit_loading(model, data_lo130, [0.5, 1, 0])
const fit_lo300 = fit_loading(model, data_lo300, [0.5, 1, 0])
const fit_lo600 = fit_loading(model, data_lo600, [0.5, 1, 0])
const fit_lo1200 = fit_loading(model, data_lo1200, [0.5, 1, 0])
const fit_lo2400 = fit_loading(model, data_lo2400, [0.5, 1, 0])
const fit_lo4800 = fit_loading(model, data_lo4800, [0.5, 1, 0])
const fit_lo9600 = fit_loading(model, data_lo9600, [0.5, 1, 0])
const fit_lo14400 = fit_loading(model, data_lo14400, [0.5, 1, 0])
const fit_lo19200 = fit_loading(model, data_lo19200, [0.5, 1, 0])

const lo_ramsey_t = [0.13, 0.3, 0.6, 1.2, 2.4, 4.8, 9.6, 14.4, 19.2]
const lo_ramsey_v = [fit_lo130.param[2],
                     fit_lo300.param[2],
                     fit_lo600.param[2],
                     fit_lo1200.param[2],
                     fit_lo2400.param[2],
                     fit_lo4800.param[2],
                     fit_lo9600.param[2],
                     fit_lo14400.param[2],
                     fit_lo19200.param[2]]
const lo_ramsey_a = [fit_lo130.unc[2],
                     fit_lo300.unc[2],
                     fit_lo600.unc[2],
                     fit_lo1200.unc[2],
                     fit_lo2400.unc[2],
                     fit_lo4800.unc[2],
                     fit_lo9600.unc[2],
                     fit_lo14400.unc[2],
                     fit_lo19200.unc[2]]

const fit_hi130 = fit_loading(model, data_hi130, [0.5, 1, 0])
const fit_hi300 = fit_loading(model, data_hi300, [0.5, 1, 0])
const fit_hi600 = fit_loading(model, data_hi600, [0.5, 1, 0])
const fit_hi1200 = fit_loading(model, data_hi1200, [0.5, 1, 0])
const fit_hi2400 = fit_loading(model, data_hi2400, [0.5, 1, 0])
const fit_hi9600 = fit_loading(model, data_hi9600, [0.5, 1, 0])

const hi_ramsey_t = [0.13, 0.3, 0.6, 1.2, 2.4]
const hi_ramsey_v = [fit_hi130.param[2],
                     fit_hi300.param[2],
                     fit_hi600.param[2],
                     fit_hi1200.param[2],
                     fit_hi2400.param[2]]
const hi_ramsey_a = [fit_hi130.unc[2],
                     fit_hi300.unc[2],
                     fit_hi600.unc[2],
                     fit_hi1200.unc[2],
                     fit_hi2400.unc[2]]

function combine_all(data)
    params, value, unc = NaCsData.get_values(NaCsData.CountData([0.0], sum(data.values.counts, dims=1)))
    return Unc(value[1, 1], unc[1, 1])
end

const aux_comb_hi130 = combine_all(aux_hi130)
const aux_comb_hi300 = combine_all(aux_hi300)
const aux_comb_hi600 = combine_all(aux_hi600)
const aux_comb_hi1200 = combine_all(aux_hi1200)
const aux_comb_hi2400 = combine_all(aux_hi2400)
const aux_comb_hi9600 = combine_all(aux_hi9600)

const hi_aux_t = [0.13, 0.3, 0.6, 1.2, 2.4]
const hi_aux_v = [aux_comb_hi130.a,
                  aux_comb_hi300.a,
                  aux_comb_hi600.a,
                  aux_comb_hi1200.a,
                  aux_comb_hi2400.a]
const hi_aux_a = [aux_comb_hi130.s,
                  aux_comb_hi300.s,
                  aux_comb_hi600.s,
                  aux_comb_hi1200.s,
                  aux_comb_hi2400.s]

function model_exp(x, p)
    return p[1] .* exp.(.-x ./ p[2])
end
function model_exp_off(x, p)
    return p[1] .* exp.(.-x ./ p[2]) .+ p[3]
end

const fit_lo_ramsey = fit_data(model_exp, lo_ramsey_t, lo_ramsey_v,
                               [0.9, 8], plot_hi=22)
const fit_hi_ramsey = fit_data(model_exp, hi_ramsey_t, hi_ramsey_v,
                               [0.9, 6], plot_hi=22)
const fit_hi_aux = fit_data(model_exp, hi_aux_t, hi_aux_v,
                            [0.9, 6], plot_hi=22)

@show fit_lo_ramsey.uncs
@show fit_hi_ramsey.uncs
@show fit_hi_aux.uncs

figure()
NaCsPlot.plot_loading_data(data_lo130, xscale=360, fmt="C0o",
         label="\$|0\\rangle_{\\mathrm{aux}}\$")
plot(fit_lo130.plotx .* 360, fit_lo130.ploty, color="C0")
NaCsPlot.plot_loading_data(data_hi130, xscale=360, fmt="C1o",
         label="\$|1\\rangle_{\\mathrm{aux}}\$")
plot(fit_hi130.plotx .* 360, fit_hi130.ploty, color="C1")
legend(fontsize=15, loc="upper center")
grid()
ylim([0, 1])
xlim([0, 360])
xlabel("Ramsey phase (\$^\\circ\$)")
ylabel("\$P_1\$")
NaCsPlot.maybe_save("$(prefix)_130us")

figure()
errorbar(lo_ramsey_t, lo_ramsey_v, lo_ramsey_a, color="C0", fmt="o",
         label="\$|0\\rangle_{\\mathrm{aux}}\$")
plot(fit_lo_ramsey.plotx, fit_lo_ramsey.ploty, color="C0")
errorbar(hi_ramsey_t, hi_ramsey_v, hi_ramsey_a, color="C1", fmt="o",
         label="\$|1\\rangle_{\\mathrm{aux}}\$")
plot(fit_hi_ramsey.plotx, fit_hi_ramsey.ploty, color="C1")
text(9, 0.38, "\$\\tau_{|0\\rangle}=$(fit_lo_ramsey.uncs[2])\$ ms", fontsize=21, color="C0")
text(6.5, 0.09, "\$\\tau_{|1\\rangle}=$(fit_hi_ramsey.uncs[2])\$ ms", fontsize=21, color="C1")
legend(fontsize=15, loc="upper right")
grid()
xlim([0, 21])
ylim([0, 0.9])
xlabel("Detection time (ms)")
ylabel("Ramsey contrast")
NaCsPlot.maybe_save("$(prefix)_ramsey")

figure()
errorbar(hi_aux_t, hi_aux_v, hi_aux_a, color="C0", fmt="o")
plot(fit_hi_aux.plotx, fit_hi_aux.ploty, color="C0")
text(1.0, 0.34, "\$\\tau_{\\mathrm{aux}}=$(fit_hi_aux.uncs[2])\$ ms", fontsize=21, color="C0")
grid()
xlim([0, 3.5])
ylim([0, 1.0])
xlabel("Detection time (ms)")
ylabel("Aux ion state")
NaCsPlot.maybe_save("$(prefix)_aux")

NaCsPlot.maybe_show()

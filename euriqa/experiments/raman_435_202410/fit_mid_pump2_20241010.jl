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

const prefix = joinpath(@__DIR__, "imgs", "data_20241010_pump_blackman")

const inames = ["000077505-Raman435MidCircuitScan.h5"]
const datas = [NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data", iname))
               for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [(range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21),
                range(0, 1, 21))]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(param .+ 1, logical,
                                               selector, maxcnt), spec)
     for ((param_name, (param, logical)), maxcnt, spec) in zip(datas, maxcnts, specs)]

function select_nth(n)
    return function selector(logicals)
        return [1, logicals[1, n]]
    end
end

function model_ramsey(x, p)
    return p[1] .+ p[2] .* sin.(2π .* (x .- p[3]))
end
function model_rabi01(x, p)
    minv = p[1]
    maxv = 1 - p[2]
    mid = (minv + maxv) / 2
    amp = (maxv - minv) / 2
    return mid .- amp .* cos.(2π .* x)
end
function model_rabi1(x, p)
    return p[1] .+ p[2] .* (cos.(2π .* x) .+ 1) / 2
end
function model_lin(x, p)
    return p[1] .+ x .* p[2]
end

const datas_0 = select_datas(datas, select_nth(2), maxcnts, specs)
const datas_1 = select_datas(datas, select_nth(3), maxcnts, specs)

const data_rabi_0 = [datas_0[1][1], datas_0[1][3], datas_0[1][5], datas_0[1][7]]
const data_rabi_1 = [datas_0[1][2], datas_0[1][4], datas_0[1][6], datas_0[1][8]]

const fit_data_rabi_0 = [fit_loading(model_rabi01, data, [0.0, 0.0])
                         for data in data_rabi_0]
const fit_data_rabi_1 = [fit_loading(model_rabi01, data, [0.0, 0.0])
                         for data in data_rabi_1]

const data_ramsey_0 = [datas_0[1][9], datas_0[1][11], datas_0[1][13], datas_0[1][15]]
const data_ramsey_1 = [datas_0[1][10], datas_0[1][12], datas_0[1][14], datas_0[1][16]]

const fit_data_ramsey_0 = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.0])
                           for data in data_ramsey_0]
const fit_data_ramsey_1 = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.0])
                           for data in data_ramsey_1]

const aux_rabi_0 = [datas_1[1][17], datas_1[1][19], datas_1[1][21]]
const aux_rabi_1 = [datas_1[1][18], datas_1[1][20]]

const fit_aux_rabi_0 = [fit_loading(model_rabi1, data, [0.0, 1.0])
                        for data in aux_rabi_0]
const fit_aux_rabi_1 = [fit_loading(model_rabi1, data, [0.0, 1.0])
                        for data in aux_rabi_1]

const pump_cycles = [0, 7, 14, 21]

const data_rabi_err0_0 = [fit.param[1] for fit in fit_data_rabi_0]
const data_rabi_err0_unc_0 = [fit.unc[1] for fit in fit_data_rabi_0]
const data_rabi_err0_1 = [fit.param[1] for fit in fit_data_rabi_1]
const data_rabi_err0_unc_1 = [fit.unc[1] for fit in fit_data_rabi_1]
const data_rabi_err1_0 = [fit.param[2] for fit in fit_data_rabi_0]
const data_rabi_err1_unc_0 = [fit.unc[2] for fit in fit_data_rabi_0]
const data_rabi_err1_1 = [fit.param[2] for fit in fit_data_rabi_1]
const data_rabi_err1_unc_1 = [fit.unc[2] for fit in fit_data_rabi_1]

const data_rabi_err0_cycle_0 = fit_data(model_lin, pump_cycles, data_rabi_err0_0,
                                        data_rabi_err0_unc_0, [0.0, 0.01])
const data_rabi_err0_cycle_1 = fit_data(model_lin, pump_cycles, data_rabi_err0_1,
                                        data_rabi_err0_unc_1, [0.0, 0.01])
const data_rabi_err1_cycle_0 = fit_data(model_lin, pump_cycles, data_rabi_err1_0,
                                        data_rabi_err1_unc_0, [0.0, 0.01])
const data_rabi_err1_cycle_1 = fit_data(model_lin, pump_cycles, data_rabi_err1_1,
                                        data_rabi_err1_unc_1, [0.0, 0.01])

const data_ramsey_err_0 = [1 - 2 * abs(fit.param[2]) for fit in fit_data_ramsey_0]
const data_ramsey_err_unc_0 = [2 * fit.unc[2] for fit in fit_data_ramsey_0]
const data_ramsey_err_1 = [1 - 2 * abs(fit.param[2]) for fit in fit_data_ramsey_1]
const data_ramsey_err_unc_1 = [2 * fit.unc[2] for fit in fit_data_ramsey_1]

const data_ramsey_err_cycle_0 = fit_data(model_lin, pump_cycles, data_ramsey_err_0,
                                         data_ramsey_err_unc_0, [0.0, 0.01])
const data_ramsey_err_cycle_1 = fit_data(model_lin, pump_cycles, data_ramsey_err_1,
                                         data_ramsey_err_unc_1, [0.0, 0.01])

const aux_rabi_err_0 = [fit.param[1] for fit in fit_aux_rabi_0]
const aux_rabi_err_unc_0 = [fit.unc[1] for fit in fit_aux_rabi_0]
const aux_rabi_err_1 = [fit.param[1] for fit in fit_aux_rabi_1]
const aux_rabi_err_unc_1 = [fit.unc[1] for fit in fit_aux_rabi_1]

figure(figsize=[6.4 * 2, 4.8 * 2])
for i in 1:4
    subplot(2, 2, i)
    fit0 = fit_data_rabi_0[i]
    plot(fit0.plotx .* 2π, fit0.ploty, color="C0")
    NaCsPlot.plot_loading_data(data_rabi_0[i], xscale=2π, fmt="C0o",
                               label="\$|0\\rangle\$")
    fit1 = fit_data_rabi_1[i]
    plot(fit1.plotx .* 2π, fit1.ploty, color="C1")
    NaCsPlot.plot_loading_data(data_rabi_1[i], xscale=2π, fmt="C1o",
                               label="\$|1\\rangle\$")
    title("$(pump_cycles[i]) cycles")
    grid()
    legend(fontsize=15)
    xlim([0, 2π])
    ylim([0, 1])
    xlabel("\$\\phi\$ (rad)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_data_rabi")

figure(figsize=[6.4 * 2, 4.8 * 2])
for i in 1:4
    subplot(2, 2, i)
    fit0 = fit_data_ramsey_0[i]
    plot(fit0.plotx .* 2π, fit0.ploty, color="C0")
    NaCsPlot.plot_loading_data(data_ramsey_0[i], xscale=2π, fmt="C0o",
                               label="\$|0\\rangle\$")
    fit1 = fit_data_ramsey_1[i]
    plot(fit1.plotx .* 2π, fit1.ploty, color="C1")
    NaCsPlot.plot_loading_data(data_ramsey_1[i], xscale=2π, fmt="C1o",
                               label="\$|1\\rangle\$")
    title("$(pump_cycles[i]) cycles")
    grid()
    legend(fontsize=15)
    xlim([0, 2π])
    ylim([0, 1])
    xlabel("\$\\phi\$ (rad)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_data_ramsey")

figure(figsize=[6.4 * 2, 4.8 * 2])
for i in 1:3
    subplot(2, 2, i)
    fit0 = fit_aux_rabi_0[i]
    plot(fit0.plotx .* 2π, fit0.ploty, color="C0")
    NaCsPlot.plot_loading_data(aux_rabi_0[i], xscale=2π, fmt="C0o",
                               label="\$|0\\rangle\$")
    if i < 3
        fit1 = fit_aux_rabi_1[i]
        plot(fit1.plotx .* 2π, fit1.ploty, color="C1")
        NaCsPlot.plot_loading_data(aux_rabi_1[i], xscale=2π, fmt="C1o",
                                   label="\$|1\\rangle\$")
    end
    title("$(pump_cycles[i]) cycles")
    grid()
    legend(fontsize=15)
    xlim([0, 2π])
    ylim([0, 1])
    xlabel("\$\\phi\$ (rad)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_aux_rabi")

figure()
plot(data_rabi_err1_cycle_0.plotx, data_rabi_err1_cycle_0.ploty, color="C0")
errorbar(pump_cycles, data_rabi_err1_0, data_rabi_err1_unc_0, fmt="C0o",
         label="\$1\\rightarrow0,|0\\rangle\$")
plot(data_rabi_err1_cycle_1.plotx, data_rabi_err1_cycle_1.ploty, color="C1")
errorbar(pump_cycles, data_rabi_err1_1, data_rabi_err1_unc_1, fmt="C1o",
         label="\$1\\rightarrow0,|1\\rangle\$")
plot(data_rabi_err0_cycle_0.plotx, data_rabi_err0_cycle_0.ploty, color="C2")
errorbar(pump_cycles, data_rabi_err0_0, data_rabi_err0_unc_0, fmt="C2o",
         label="\$0\\rightarrow1,|0\\rangle\$")
plot(data_rabi_err0_cycle_1.plotx, data_rabi_err0_cycle_1.ploty, color="C3")
errorbar(pump_cycles, data_rabi_err0_1, data_rabi_err0_unc_1, fmt="C3o",
         label="\$0\\rightarrow1,|1\\rangle\$")
text(9, 0.12, "\$E_{1\\rightarrow0,|0\\rangle}=$(data_rabi_err1_cycle_0.uncs[2] * 100) \\%\$",
     color="C0", fontsize=15)
text(9, 0.106, "\$E_{1\\rightarrow0,|1\\rangle}=$(data_rabi_err1_cycle_1.uncs[2] * 100) \\%\$",
     color="C1", fontsize=15)
text(9, 0.092, "\$E_{0\\rightarrow1,|0\\rangle}=$(data_rabi_err0_cycle_0.uncs[2] * 100) \\%\$",
     color="C2", fontsize=15)
text(9, 0.078, "\$E_{0\\rightarrow1,|1\\rangle}=$(data_rabi_err0_cycle_1.uncs[2] * 100) \\%\$",
     color="C3", fontsize=15)
xlim([0, 23])
ylim([0, 0.13])
legend(fontsize=15)
grid()
NaCsPlot.maybe_save("$(prefix)_data_rabi_error")

figure()
plot(data_ramsey_err_cycle_0.plotx, data_ramsey_err_cycle_0.ploty, color="C0")
errorbar(pump_cycles, data_ramsey_err_0, data_ramsey_err_unc_0, fmt="C0o",
         label="\$|0\\rangle\$")
plot(data_ramsey_err_cycle_1.plotx, data_ramsey_err_cycle_1.ploty, color="C1")
errorbar(pump_cycles, data_ramsey_err_1, data_ramsey_err_unc_1, fmt="C1o",
         label="\$|1\\rangle\$")
text(10, 0.01, "\$E_{|0\\rangle}=$(data_ramsey_err_cycle_0.uncs[2] * 100) \\%\$",
     color="C0", fontsize=16)
text(10, 0.05, "\$E_{|1\\rangle}=$(data_ramsey_err_cycle_1.uncs[2] * 100) \\%\$",
     color="C1", fontsize=16)
xlim([0, 23])
ylim([0, 0.30])
legend(fontsize=15)
grid()
NaCsPlot.maybe_save("$(prefix)_data_ramsey_error")

figure()
errorbar(pump_cycles[1:3], aux_rabi_err_0, aux_rabi_err_unc_0, fmt="C0o-",
         label="\$|0\\rangle\$")
errorbar(pump_cycles[1:2], aux_rabi_err_1, aux_rabi_err_unc_1, fmt="C1o-",
         label="\$|1\\rangle\$")
xlim([0, 23])
ylim([0, 0.06])
grid()
NaCsPlot.maybe_save("$(prefix)_aux_error")

NaCsPlot.maybe_show()

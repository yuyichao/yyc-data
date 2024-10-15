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

const prefix = joinpath(@__DIR__, "imgs", "data_20241010_pump_77469")

const inames = ["000077469-Raman435MidCircuitScan.h5"]
const datas = [NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data", iname))
               for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [(range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11))]

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
function model_rabi(x, p)
    return p[1] .+ p[2] .* (cos.(2π .* x) .+ 1) / 2
end

const datas_0 = select_datas(datas, select_nth(2), maxcnts, specs)
const datas_1 = select_datas(datas, select_nth(3), maxcnts, specs)

const data_ramsey_0 = [datas_0[1][1], datas_0[1][3], datas_0[1][5], datas_0[1][7]]
const data_ramsey_1 = [datas_0[1][2], datas_0[1][4], datas_0[1][6], datas_0[1][8]]

const fit_data_ramsey_0 = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.0])
                           for data in data_ramsey_0]
const fit_data_ramsey_1 = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.0])
                           for data in data_ramsey_1]

const aux_rabi_0 = [datas_1[1][9], datas_1[1][11], datas_1[1][13], datas_1[1][15]]
const aux_rabi_1 = [datas_1[1][10], datas_1[1][12], datas_1[1][14], datas_1[1][16]]

const fit_aux_rabi_0 = [fit_loading(model_rabi, data, [0.0, 1.0])
                        for data in aux_rabi_0]
const fit_aux_rabi_1 = [fit_loading(model_rabi, data, [0.0, 1.0])
                        for data in aux_rabi_1]

const pump_cycles = [0, 7, 14, 21]

const data_ramsey_err_0 = [1 - 2 * abs(fit.param[2]) for fit in fit_data_ramsey_0]
const data_ramsey_err_unc_0 = [2 * fit.unc[2] for fit in fit_data_ramsey_0]
const data_ramsey_err_1 = [1 - 2 * abs(fit.param[2]) for fit in fit_data_ramsey_1]
const data_ramsey_err_unc_1 = [2 * fit.unc[2] for fit in fit_data_ramsey_1]

const aux_rabi_err_0 = [fit.param[1] for fit in fit_aux_rabi_0]
const aux_rabi_err_unc_0 = [fit.unc[1] for fit in fit_aux_rabi_0]
const aux_rabi_err_1 = [fit.param[1] for fit in fit_aux_rabi_1]
const aux_rabi_err_unc_1 = [fit.unc[1] for fit in fit_aux_rabi_1]

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
    legend()
    xlim([0, 2π])
    ylim([0, 1])
    xlabel("\$\\phi\$ (rad)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_data_ramsey")

figure(figsize=[6.4 * 2, 4.8 * 2])
for i in 1:4
    subplot(2, 2, i)
    fit0 = fit_aux_rabi_0[i]
    plot(fit0.plotx .* 2π, fit0.ploty, color="C0")
    NaCsPlot.plot_loading_data(aux_rabi_0[i], xscale=2π, fmt="C0o",
                               label="\$|0\\rangle\$")
    fit1 = fit_aux_rabi_1[i]
    plot(fit1.plotx .* 2π, fit1.ploty, color="C1")
    NaCsPlot.plot_loading_data(aux_rabi_1[i], xscale=2π, fmt="C1o",
                               label="\$|1\\rangle\$")
    title("$(pump_cycles[i]) cycles")
    grid()
    legend()
    xlim([0, 2π])
    ylim([0, 1])
    xlabel("\$\\phi\$ (rad)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_aux_rabi")

figure()
errorbar(pump_cycles, data_ramsey_err_0, data_ramsey_err_unc_0, fmt="C0o-")
errorbar(pump_cycles, data_ramsey_err_1, data_ramsey_err_unc_1, fmt="C2o-")
errorbar(pump_cycles, aux_rabi_err_0, aux_rabi_err_unc_0, fmt="C1o-")
errorbar(pump_cycles, aux_rabi_err_1, aux_rabi_err_unc_1, fmt="C3o-")
xlim([0, 23])
ylim([0, 0.25])
grid()
NaCsPlot.maybe_save("$(prefix)_error")

NaCsPlot.maybe_show()

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

const prefix = joinpath(@__DIR__, "imgs", "data_20241017_pump_echo_slow")

const inames = ["000079495-Raman435MidCircuitScan.h5",
                "000079497-Raman435MidCircuitScan.h5"]
const datas = [NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data", iname))
               for iname in inames]
const maxcnts = [typemax(Int), typemax(Int)]
const specs = [(([0], [0], [8], [8], [16], [16], [24], [24]),
                ([0], [0], [8], [8], [16], [16], [24], [24]),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11),
                range(0, 1, 11)),
               (range(0, 1, 11),
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

const pump_cycles = [0, 8, 16, 24]

const data_z_err0_aux0 = [datas_0[1][1][1]; datas_0[1][1][3];
                          datas_0[1][1][5]; datas_0[1][1][7]]
const data_z_err0_aux1 = [datas_0[1][1][2]; datas_0[1][1][4];
                          datas_0[1][1][6]; datas_0[1][1][8]]
const data_z_err1_aux0 = [datas_0[1][2][1]; datas_0[1][2][3];
                          datas_0[1][2][5]; datas_0[1][2][7]]
const data_z_err1_aux1 = [datas_0[1][2][2]; datas_0[1][2][4];
                          datas_0[1][2][6]; datas_0[1][2][8]]

const fit_data_z_err0_aux0 = fit_loading(model_lin, data_z_err0_aux0, [0.0, 0.01])
const fit_data_z_err0_aux1 = fit_loading(model_lin, data_z_err0_aux1, [0.0, 0.01])
const fit_data_z_err1_aux0 = fit_loading(model_lin, data_z_err1_aux0, [0.0, 0.01])
const fit_data_z_err1_aux1 = fit_loading(model_lin, data_z_err1_aux1, [0.0, 0.01])

const data_ramsey_aux0 = [datas_0[1][3], datas_0[1][5], datas_0[1][7], datas_0[1][9]]
const data_ramsey_aux1 = [datas_0[1][4], datas_0[1][6], datas_0[1][8], datas_0[1][10]]

const fit_data_ramsey_aux0 = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.0])
                              for data in data_ramsey_aux0]
const fit_data_ramsey_aux1 = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.0])
                              for data in data_ramsey_aux1]

const data_x_err_aux0 = [1 - 2 * abs(fit.param[2]) for fit in fit_data_ramsey_aux0]
const data_x_err_unc_aux0 = [2 * fit.unc[2] for fit in fit_data_ramsey_aux0]
const data_x_err_aux1 = [1 - 2 * abs(fit.param[2]) for fit in fit_data_ramsey_aux1]
const data_x_err_unc_aux1 = [2 * fit.unc[2] for fit in fit_data_ramsey_aux1]

const fit_data_x_err_aux0 = fit_data(model_lin, pump_cycles,
                                       data_x_err_aux0,
                                       data_x_err_unc_aux0, [0.0, 0.01])
const fit_data_x_err_aux1 = fit_data(model_lin, pump_cycles,
                                       data_x_err_aux1,
                                       data_x_err_unc_aux1, [0.0, 0.01])

const aux_rabi_0 = [datas_1[2][1], datas_1[2][3], datas_1[2][5], datas_1[2][7]]
const aux_rabi_1 = [datas_1[2][2], datas_1[2][4], datas_1[2][6], datas_1[2][8]]

const fit_aux_rabi_0 = [fit_loading(model_rabi1, data, [0.0, 1.0])
                        for data in aux_rabi_0]
const fit_aux_rabi_1 = [fit_loading(model_rabi1, data, [0.0, 1.0])
                        for data in aux_rabi_1]

const aux_rabi_err_0 = [fit.param[1] for fit in fit_aux_rabi_0]
const aux_rabi_err_unc_0 = [fit.unc[1] for fit in fit_aux_rabi_0]
const aux_rabi_err_1 = [fit.param[1] for fit in fit_aux_rabi_1]
const aux_rabi_err_unc_1 = [fit.unc[1] for fit in fit_aux_rabi_1]

figure(figsize=[6.4 * 2, 4.8 * 2])
for i in 1:4
    subplot(2, 2, i)
    fit0 = fit_data_ramsey_aux0[i]
    plot(fit0.plotx .* 2π, fit0.ploty, color="C0")
    NaCsPlot.plot_loading_data(data_ramsey_aux0[i], xscale=2π, fmt="C0o",
                               label="\$|0\\rangle\$")
    fit1 = fit_data_ramsey_aux1[i]
    plot(fit1.plotx .* 2π, fit1.ploty, color="C1")
    NaCsPlot.plot_loading_data(data_ramsey_aux1[i], xscale=2π, fmt="C1o",
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
    legend(fontsize=15)
    xlim([0, 2π])
    ylim([0, 1])
    xlabel("\$\\phi\$ (rad)")
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_aux_rabi")

figure()
plot(fit_data_z_err1_aux0.plotx, fit_data_z_err1_aux0.ploty, color="C0")
NaCsPlot.plot_loading_data(data_z_err1_aux0, fmt="C0o",
         label="\$1\\rightarrow0,|0\\rangle\$")
plot(fit_data_z_err1_aux1.plotx, fit_data_z_err1_aux1.ploty, color="C1")
NaCsPlot.plot_loading_data(data_z_err1_aux1, fmt="C1o",
         label="\$1\\rightarrow0,|1\\rangle\$")
plot(fit_data_z_err0_aux0.plotx, fit_data_z_err0_aux0.ploty, color="C2")
NaCsPlot.plot_loading_data(data_z_err0_aux0, fmt="C2o",
         label="\$0\\rightarrow1,|0\\rangle\$")
plot(fit_data_z_err0_aux1.plotx, fit_data_z_err0_aux1.ploty, color="C3")
NaCsPlot.plot_loading_data(data_z_err0_aux1, fmt="C3o",
         label="\$0\\rightarrow1,|1\\rangle\$")
text(10, 0.036, "\$E_{1\\rightarrow0,|0\\rangle}=$(fit_data_z_err1_aux0.uncs[2] * 100) \\%\$/cycle",
     color="C0", fontsize=15)
text(10, 0.0318, "\$E_{1\\rightarrow0,|1\\rangle}=$(fit_data_z_err1_aux1.uncs[2] * 100) \\%\$/cycle",
     color="C1", fontsize=15)
text(10, 0.0276, "\$E_{0\\rightarrow1,|0\\rangle}=$(fit_data_z_err0_aux0.uncs[2] * 100) \\%\$/cycle",
     color="C2", fontsize=15)
text(10, 0.0234, "\$E_{0\\rightarrow1,|1\\rangle}=$(fit_data_z_err0_aux1.uncs[2] * 100) \\%\$/cycle",
     color="C3", fontsize=15)
xlim([0, 25])
ylim([0, 0.039])
legend(fontsize=15)
grid()
NaCsPlot.maybe_save("$(prefix)_data_z_error")

figure()
plot(fit_data_x_err_aux0.plotx, fit_data_x_err_aux0.ploty .* 100, color="C0")
errorbar(pump_cycles, data_x_err_aux0 .* 100, data_x_err_unc_aux0 .* 100, fmt="C0o",
         label="\$|0\\rangle\$")
plot(fit_data_x_err_aux1.plotx, fit_data_x_err_aux1.ploty .* 100, color="C1")
errorbar(pump_cycles, data_x_err_aux1 .* 100, data_x_err_unc_aux1 .* 100, fmt="C1o",
         label="\$|1\\rangle\$")
text(10, 0.2, "\$E_{|0\\rangle}=$(fit_data_x_err_aux0.uncs[2] * 100) \\%\$/cycle",
     color="C0", fontsize=16)
text(10, 1.2, "\$E_{|1\\rangle}=$(fit_data_x_err_aux1.uncs[2] * 100) \\%\$/cycle",
     color="C1", fontsize=16)
xlim([0, 25])
ylim([0, 8])
xlabel("Pumping cycles")
ylabel("Error (%)")
legend(fontsize=15)
grid()
NaCsPlot.maybe_save("$(prefix)_data_x_error")

figure()
errorbar(pump_cycles, aux_rabi_err_0 .* 100, aux_rabi_err_unc_0 .* 100, fmt="C0o-",
         label="\$|0\\rangle\$")
errorbar(pump_cycles, aux_rabi_err_1 .* 100, aux_rabi_err_unc_1 .* 100, fmt="C1o-",
         label="\$|1\\rangle\$")
xlim([0, 25])
ylim([0, 6])
xlabel("Pumping cycles")
ylabel("Error (%)")
grid()
NaCsPlot.maybe_save("$(prefix)_aux_error")

NaCsPlot.maybe_show()

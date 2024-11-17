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

const prefix = joinpath(@__DIR__, "imgs", "data_20240930_detect")

const inames = ["000074901-Raman435MidCircuitScan.h5",
                "000074914-Raman435MidCircuitScan.h5",
                "000074948-Raman435MidCircuitScan.h5",
                "000075679-Raman435MidCircuitScan.h5",
                "000075680-Raman435MidCircuitScan.h5"]
const datas = [NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data", iname),
                                                index_param=true)
               for iname in inames]
const maxcnts = [typemax(Int), typemax(Int), typemax(Int), typemax(Int), typemax(Int)]
const specs = [range(0, 1, 21),
               range(0, 1, 21),
               range(0, 1, 21),
               range(0, 1, 21),
               range(0, 1, 21)]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(param, logical,
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

const fit_data_ramsey = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.0])
                         for data in datas_0]
const fit_aux_rabi = fit_loading(model_rabi01, datas_1[3], [0.0, 0.0])

figure()
plot(fit_data_ramsey[1].plotx .* 2π, fit_data_ramsey[1].ploty, color="C0")
NaCsPlot.plot_loading_data(datas_0[1], xscale=2π, fmt="C0o",
                           label="\$|0\\rangle\$")
plot(fit_data_ramsey[2].plotx .* 2π, fit_data_ramsey[2].ploty, color="C1")
NaCsPlot.plot_loading_data(datas_0[2], xscale=2π, fmt="C1o",
                           label="\$|1\\rangle\$")
text(1, 0.22, "\$E_{|0\\rangle}=$((1 - fit_data_ramsey[1].uncs[2] * 2) * 100) \\%\$",
     color="C0", fontsize=15)
text(1, 0.12, "\$E_{|1\\rangle}=$((1 - fit_data_ramsey[2].uncs[2] * 2) * 100) \\%\$",
     color="C1", fontsize=15)
xlim([0, 2π])
ylim([0, 1])
xlabel("\$\\phi\$ (rad)")
grid()
legend(fontsize=15)
NaCsPlot.maybe_save("$(prefix)_data_ramsey_square")

figure()
plot(fit_data_ramsey[4].plotx .* 2π, fit_data_ramsey[4].ploty, color="C0")
NaCsPlot.plot_loading_data(datas_0[4], xscale=2π, fmt="C0o",
                           label="\$|0\\rangle\$")
plot(fit_data_ramsey[5].plotx .* 2π, fit_data_ramsey[5].ploty, color="C1")
NaCsPlot.plot_loading_data(datas_0[5], xscale=2π, fmt="C1o",
                           label="\$|1\\rangle\$")
text(1, 0.22, "\$E_{|0\\rangle}=$((1 - fit_data_ramsey[4].uncs[2] * 2) * 100) \\%\$",
     color="C0", fontsize=15)
text(1, 0.12, "\$E_{|1\\rangle}=$((1 - fit_data_ramsey[5].uncs[2] * 2) * 100) \\%\$",
     color="C1", fontsize=15)
xlim([0, 2π])
ylim([0, 1])
xlabel("\$\\phi\$ (rad)")
grid()
legend(fontsize=15)
NaCsPlot.maybe_save("$(prefix)_data_ramsey_half_sk1")

println(fit_data_ramsey[1].uncs[2] * 2)
println(fit_data_ramsey[2].uncs[2] * 2)
println(fit_data_ramsey[4].uncs[2] * 2)
println(fit_data_ramsey[5].uncs[2] * 2)

println(fit_aux_rabi.uncs)

figure()
plot(fit_aux_rabi.plotx .* 2π, fit_aux_rabi.ploty, color="C0")
NaCsPlot.plot_loading_data(datas_1[3], xscale=2π, fmt="C0o")
text(2.05, 0.22, "\$E_{|0\\rangle}=$(fit_aux_rabi.uncs[1] * 100) \\%\$",
     color="C0", fontsize=15)
text(2.05, 0.12, "\$E_{|1\\rangle}=$(fit_aux_rabi.uncs[2] * 100) \\%\$",
     color="C0", fontsize=15)
xlim([0, 2π])
ylim([0, 1])
xlabel("\$\\theta\$ (rad)")
grid()
NaCsPlot.maybe_save("$(prefix)_aux_rabi")

NaCsPlot.maybe_show()

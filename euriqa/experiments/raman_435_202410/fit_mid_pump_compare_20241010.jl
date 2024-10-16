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

const prefix = joinpath(@__DIR__, "imgs", "data_20241010_pump_compare")

const inames = ["000077487-Raman435MidCircuitScan.h5"]
const datas = [NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data", iname))
               for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [((range(0, 1, 11),
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
                 range(0, 1, 11)),
                (range(0, 1, 11),
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
                 range(0, 1, 11)))]

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

const data_rabi_0 = select_datas(datas, select_nth(2), maxcnts, specs)[1][1]
const fit_rabi_0 = [fit_loading(model_rabi01, data, [0.0, 0.0])
                    for data in data_rabi_0]
const data_ramsey_0 = select_datas(datas, select_nth(2), maxcnts, specs)[1][2]
const fit_ramsey_0 = [fit_loading(model_ramsey, data, [0.5, 0.5, 0.1])
                      for data in data_ramsey_0]

const amps = [abs(fit.param[2]) for fit in fit_ramsey_0]
const amp_errs = [fit.unc[2] for fit in fit_ramsey_0]

function real_model(x, p)
    val = p[1] / 2
    idx = Int(x)
    if idx & 1 != 0
        val *= p[2]
    end
    if idx & 2 != 0
        val *= p[3]
    end
    if idx & 4 != 0
        val *= p[4]
    end
    if idx & 8 != 0
        val *= p[5]
    end
    return val
end

function model(idx, p)
    return real_model.(idx, Ref(p))
end

const fit_err = fit_data(model, 0:15, amps, amp_errs, [1.0, 1.0, 1.0, 1.0, 1.0],
                         plotx=0:15)
@show fit_err
const plot_names = ["none",
                    "-1",
                    "1",
                    "-1,1",
                    "0",
                    "-1,0",
                    "0,1",
                    "-1,0,1",
                    "935",
                    "935,-1",
                    "935,1",
                    "935,-1,1",
                    "935,0",
                    "935,-1,0",
                    "935,0,1",
                    "935,-1,0,1"]

figure(figsize=[6.4 * 4, 4.8 * 4])
for i in 1:16
    data = data_rabi_0[i]
    fit = fit_rabi_0[i]
    subplot(4, 4, i)
    plot(fit.plotx .* 2π, fit.ploty, color="C0")
    NaCsPlot.plot_loading_data(data, xscale=2π, fmt="C0o")
    text(2.25, 0.52, "$(100 - fit.uncs[2] * 100) %", color="C0")
    text(2.25, 0.10, "$(fit.uncs[1] * 100) %", color="C0")
    title(plot_names[i])
    grid()
    xlim([0, 2π])
    ylim([0, 1])
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_data_rabi")

figure(figsize=[6.4 * 4, 4.8 * 4])
for i in 1:16
    data = data_ramsey_0[i]
    fit = fit_ramsey_0[i]
    subplot(4, 4, i)
    plot(fit.plotx .* 2π, fit.ploty, color="C0")
    NaCsPlot.plot_loading_data(data, xscale=2π, fmt="C0o")
    text(2, 0.1, "$(abs(fit.uncs[2]) * 200) %", color="C0")
    title(plot_names[i])
    grid()
    xlim([0, 2π])
    ylim([0, 1])
end
tight_layout()
NaCsPlot.maybe_save("$(prefix)_data_ramsey")

NaCsPlot.maybe_show()

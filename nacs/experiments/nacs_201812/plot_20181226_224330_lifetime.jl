#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181226_224330.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [OrderedDict(
    :hot=>[0, 50, 100, 200, 500, 1000, 2000, 5000, 10000],
    :cold=>[0, 50, 100, 200, 500, 1000, 2000, 5000, 10000],
)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function fit_survival(model, data, p0; plotx=nothing, plot_lo=nothing, plot_hi=nothing,
                      use_unc=true, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        if plot_lo === nothing
            plot_lo = mid - span * plot_scale / 2
            if plot_lo * lo <= 0
                plot_lo = 0
            end
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
            if plot_hi * hi <= 0
                plot_hi = 0
            end
        end
        plotx = linspace(plot_lo, plot_hi, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], uncs[:, 2].^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)

model_exp(x, p) = p[1] .* exp.(x ./ -p[2])

data_na_cold = datas_na[1][:cold]
data_na_hot = datas_na[1][:hot]
data_cs_cold = datas_cs[1][:cold]
data_cs_hot = datas_cs[1][:hot]

fit_na_cold = fit_survival(model_exp, data_na_cold, [0.9, 500])
fit_na_hot = fit_survival(model_exp, data_na_hot, [0.9, 500])
fit_cs_cold = fit_survival(model_exp, data_cs_cold, [0.9, 5000])
fit_cs_hot = fit_survival(model_exp, data_cs_hot, [0.9, 5000])

const prefix = joinpath(@__DIR__, "imgs", "data_20181226_224330_lifetime")

figure()
NaCsPlot.plot_survival_data(data_na_cold, fmt="C0.", label="Cold")
plot(fit_na_cold.plotx, fit_na_cold.ploty)
NaCsPlot.plot_survival_data(data_na_hot, fmt="C1.", label="Hot")
plot(fit_na_hot.plotx, fit_na_hot.ploty)
text(1500, 0.5, "\$\\tau=$(fit_na_cold.uncs[2])\$ms", color="C0")
text(1500, 0.35, "\$\\tau=$(fit_na_hot.uncs[2])\$ms", color="C1")
legend()
grid()
ylim([0, 1])
xlim([0, 3000])
title("Na Lifetime")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs_cold, fmt="C0.", label="Cold")
plot(fit_cs_cold.plotx, fit_cs_cold.ploty)
NaCsPlot.plot_survival_data(data_cs_hot, fmt="C1.", label="Hot")
plot(fit_cs_hot.plotx, fit_cs_hot.ploty)
text(5500, 0.5, "\$\\tau=$(fit_cs_cold.uncs[2] / 1000)\$s", color="C0")
text(5500, 0.35, "\$\\tau=$(fit_cs_hot.uncs[2] / 1000)\$s", color="C1")
legend()
grid()
ylim([0, 1])
xlim([0, 11000])
title("Cs Lifetime")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()

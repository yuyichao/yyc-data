#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200105_123902.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [([0, 0.1, 0.2, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3.5, 4.5, 5, 7, 10],
                368.1571 .+ [-80; -20:4:20; 80] .* 1e-3)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

fit_data(model, x, y, p0; kws...) =
    fit_data(model, x, y, nothing, p0; kws...)

function fit_data(model, params, ratios, uncs, p0;
                  plotx=nothing, plot_lo=nothing, plot_hi=nothing, plot_scale=1.1)
    use_unc = uncs !== nothing
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
        fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function fit_survival(model, data, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20200105_123902_raman_time_3311")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_expoff(x, p)
    p[1] .+ p[2] .* exp.(.- x ./ p[3])
end
fit1 = fit_survival(model_expoff, datas_nacs[1][1], [0.39, 0.39, 1])
fit2 = fit_survival(model_lorentzian, datas_nacs[1][2], [0.7, 0.1, 368.1571, 0.05])

# @show fit2.uncs

figure()
NaCsPlot.plot_survival_data(datas_nacs[1][1], fmt="C0.")
plot(fit1.plotx, fit1.ploty, "C0-")
xlim([0, 10.5])
text(4, 0.58, "\$\\tau=$(fit1.uncs[3])\$ ms", color="C0")
# legend(fontsize="small", loc="lower right")
title("306460.7 GHz, 8 mw, 368.1571 MHz")
grid()
xlabel("Time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_time")

figure()
NaCsPlot.plot_survival_data(datas_nacs[1][2], fmt="C0.", label="1.25 ms")
plot(fit2.plotx, fit2.ploty, "C0-")
ylim([0.58, 0.72])
text(368.08, 0.6, "\$f=$(fit2.uncs[3])\$ MHz", color="C0", fontsize="small")
text(368.08, 0.62, "\$\\Gamma=$(fit2.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
# legend(fontsize="small", loc="lower right")
title("306460.7 GHz, 8 mw, 1.25 ms")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_det")

NaCsPlot.maybe_show()

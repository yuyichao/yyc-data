#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200208_004646.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [(308.3 .+ [-3000; -800:100:800; 2000] .* 1e-3,
                308.3 .+ [-3000; -800:100:800; 2000] .* 1e-3)]
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

const data_1 = datas_nacs[1][1]
const data_03 = datas_nacs[1][2]

const prefix = joinpath(@__DIR__, "imgs", "data_20200208_004646_raman_det")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_1 = fit_survival(model_lorentzian, data_1, [0.7, 0.2, 308.3, 0.01])
fit_03 = fit_survival(model_lorentzian, data_03, [0.7, 0.2, 307.8, 0.01])

figure()
NaCsPlot.plot_survival_data(data_1, fmt="C0.", label="ratio 1, 0.5 ms")
plot(fit_1.plotx, fit_1.ploty, "C0-")
NaCsPlot.plot_survival_data(data_03, fmt="C1.", label="raito 0.03, 1 ms")
plot(fit_03.plotx, fit_03.ploty, "C1-")
title("306491 GHz, 16 mW")
text(304.9, 0.22, "\$f_{1}=$(fit_1.uncs[3])\$ MHz", color="C0", fontsize="small")
text(304.9, 0.27, "\$\\Gamma_{1}=$(fit_1.uncs[4])\$ MHz", color="C0", fontsize="small")
text(304.9, 0.33, "\$f_{0.03}=$(fit_03.uncs[3])\$ MHz", color="C1", fontsize="small")
text(304.9, 0.38, "\$\\Gamma_{0.03}=$(fit_03.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
legend(fontsize="small", loc="lower right", handlelength=0.05)
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181213_200552.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [(0:0.75:15,
                   0:0.75:15,
                   0:0.75:15)]
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

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)

data_na = datas_na[1]

model_sin(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4])

fit_2 = fit_survival(model_sin, data_na[2], [0.4, 0.4, 0.12, -3])
fit_3 = fit_survival(model_sin, data_na[3], [0.4, 0.4, 0.32, -3])

const prefix = joinpath(@__DIR__, "imgs", "data_20181213_200552_ramsey_freq_check")

figure()
NaCsPlot.plot_survival_data(data_na[1], fmt="C0.-", label="-24 kHz")
NaCsPlot.plot_survival_data(data_na[2], fmt="C1.", label="100 kHz")
plot(fit_2.plotx, fit_2.ploty, "C1-")
NaCsPlot.plot_survival_data(data_na[3], fmt="C2.", label="300 kHz")
plot(fit_3.plotx, fit_3.ploty, "C2-")
text(7, 0.89, "\$\\nu_{100}=$(fit_2.uncs[3] * 1000)\$kHz", color="C1")
text(7, 0.78, "\$\\nu_{300}=$(fit_3.uncs[3] * 1000)\$kHz", color="C2")
legend(fontsize="small")
grid()
xlim([0, 16])
ylim([0, 1])
xlabel("Raman time (us)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190318_005439.mat",
                "data_20190318_143842.mat",
                ]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 ]
const specs = [([0],
                [0],
                [0, 0.2, 0.5, 1, 2, 4, 8, 12, 16, 24, 32],
                [0, 0.2, 0.5, 1, 2, 4, 8, 12, 16, 24, 32]),
               ([0],
                [0],
                [0, 0.2, 0.5, 1, 2, 4, 8, 12, 16, 24, 32],
                [0, 0.2, 0.5, 1, 2, 4, 8, 12, 16, 24, 32]),
               ]
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

const datas_cs = select_datas(datas, NaCsData.select_single((1, 2), (4,)), maxcnts, specs)

model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

data_488 = datas_cs[1]
data_635 = datas_cs[2]

fit_488 = fit_survival(model_exp_off, data_488[3], [0.5, 5, 0.5])
fit_635 = fit_survival(model_exp_off, data_635[3], [0.5, 5, 0.5])

# @show fit_488.uncs
# @show fit_635.uncs

const prefix = joinpath(@__DIR__, "imgs", "data_20190318_lifetime_32")

figure()
NaCsPlot.plot_survival_data(data_488[3], fmt="C0.")
plot(fit_488.plotx, fit_488.ploty, "C0-", label="488 GHz Total")
NaCsPlot.plot_survival_data(data_488[4], fmt="C1.-", label="488 GHz F = 3")
NaCsPlot.plot_survival_data(data_635[3], fmt="C2.")
plot(fit_635.plotx, fit_635.ploty, "C2-", label="635 GHz Total")
NaCsPlot.plot_survival_data(data_635[4], fmt="C3.-", label="635 GHz F = 3")
text(0.06, 0.9, "\$\\tau_{488}=$(fit_488.uncs[2]) ms\$", color="C0")
text(1.7, 0.78, "\$\\tau_{635}=$(fit_635.uncs[2]) ms\$", color="C2")
legend(fontsize="small")
grid()
ax = gca()
ax[:set_xscale]("log", nonposx="clip")
ax[:set_xticks]([0.1, 1, 10])
ax[:set_xticklabels](["0.1", "1", "10"])
ylim([0, 1])
xlim([0.05, 40])
title("Cs / Na + Cs")
xlabel("Wait time (ms)")
ylabel("Cs survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

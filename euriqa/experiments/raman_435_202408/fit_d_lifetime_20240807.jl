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

const prefix = joinpath(@__DIR__, "imgs", "data_20240807_d_lifetime")

const inames = ["000066669-DF2LifetimeScan.h5"]
const datas = [NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data", iname),
                                                index_param=true)
               for iname in inames]
const maxcnts = [typemax(Int), typemax(Int), typemax(Int)]
const specs = [range(0, 80, 17)]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(param, logical, selector, maxcnt), spec)
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
function model_exp1(x, p)
    return 1 .- p[1] .* exp.(.-x ./ p[2])
end
function model_exp(x, p)
    return p[3] .- p[1] .* exp.(.-x ./ p[2])
end

const datas_0 = select_datas(datas, select_nth(2), maxcnts, specs)

const decay_d = datas_0[1]
# const fit_decay_d = fit_loading(model_exp1, decay_d, [1.0, 10])
const fit_decay_d = fit_loading(model_exp, decay_d, [1.0, 10, 1.0])

@show fit_decay_d.uncs

figure()
plot(fit_decay_d.plotx, fit_decay_d.ploty, "C0")
NaCsPlot.plot_loading_data(decay_d, fmt="C0o")
ylim([0, 1])
xlabel("Wait time (ms)")
ylabel("S population")
grid()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

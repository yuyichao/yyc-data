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

const prefix = joinpath(@__DIR__, "imgs", "data_20241010_435_f1_spec")

const inames = ["000077388-Raman435MidCircuitScan.h5",
                "000077389-Raman435MidCircuitScan.h5",
                "000077393-Raman435MidCircuitScan.h5"]
const datas = [NaCsData.load_dax_scan_logicals1(joinpath(@__DIR__, "data", iname),
                                                index_param=true)
               for iname in inames]
const maxcnts = [typemax(Int), typemax(Int), typemax(Int)]
const specs = [range(4.1, 5.7, 101),
               range(4.89, 5.09, 41),
               range(4.1, 5.7, 201)]

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

const datas_0 = select_datas(datas, select_nth(2), maxcnts, specs)
const datas_1 = select_datas(datas, select_nth(3), maxcnts, specs)

const spec_raw = [datas_1[1]; datas_1[2]; datas_1[3]]
const spec_dressed = [datas_0[1]; datas_0[2]; datas_0[3]]

figure()
NaCsPlot.plot_loading_data(spec_raw, fmt="C0o-", label="Raw")
NaCsPlot.plot_loading_data(spec_dressed, fmt="C1o-", label="Dressed")
# xlim([-0.26, 0.06])
ylim([0, 1])
xlabel("Frequency")
ylabel("S population")
grid()
legend(fontsize=15)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

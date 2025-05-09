#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180111_123031.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec = OrderedDict(
    :x=>0:25,
    :y=>0:25,
)

const split_a = NaCsData.split_data(data_a, spec)

const data_x = split_a[:x]
const data_y = split_a[:y]

const prefix = joinpath(@__DIR__, "imgs", "data_20180111_123031_rabi_r_hot")

figure()
NaCsPlot.plot_survival_data(data_x, fmt="o-")
grid()
ylim([0, 0.7])
title("X")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y, fmt="o-")
grid()
ylim([0, 0.7])
title("Y")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

NaCsPlot.maybe_show()

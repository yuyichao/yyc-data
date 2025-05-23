#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180105_195845.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[4, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec = OrderedDict(
    :z=>linspace(-20, 100, 61),
    :x=>linspace(-180, 220, 51),
    :y=>linspace(-180, 220, 51),
)

const split_a = NaCsData.split_data(data_a, spec)

const data_z = split_a[:z]
const data_x = split_a[:x]
const data_y = split_a[:y]

const prefix = joinpath(@__DIR__, "imgs", "data_20180105_195845_hot")
const data_prefix = joinpath(@__DIR__, "data", "data_20180105_195845_hot")

figure()
NaCsPlot.plot_survival_data(data_z, fmt="o-")
grid()
ylim([0, 0.4])
title("Z spectrum")
xlabel("Frequency (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

figure()
NaCsPlot.plot_survival_data(data_x, fmt="o-")
grid()
ylim([0, 0.55])
title("X spectrum")
xlabel("Frequency (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y, fmt="o-")
grid()
ylim([0, 0.55])
title("Y spectrum")
xlabel("Frequency (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

NaCsData.dump_raw("$(data_prefix)_x.csv", data_x)
NaCsData.dump_raw("$(data_prefix)_y.csv", data_y)
NaCsData.dump_raw("$(data_prefix)_z.csv", data_z)

NaCsPlot.maybe_show()

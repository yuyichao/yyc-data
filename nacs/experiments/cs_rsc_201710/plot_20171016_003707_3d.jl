#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20171016_003707.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :axial_z=>-20:2:120,
    :radial_x=>-150:8:250,
    :radial_y=>-150:8:250,
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20171016_003707_3d")
const sorted_prefix = joinpath(@__DIR__, "sorted", "data_20171016_003707_3d")

to_sideband(f) = (i, v)->(v - f)

data_rx = NaCsData.map_params(to_sideband(46), split_a[:radial_x])
data_ry = NaCsData.map_params(to_sideband(46), split_a[:radial_y])
data_az = NaCsData.map_params(to_sideband(47), split_a[:axial_z])

NaCsData.dump_raw("$(sorted_prefix)_rx.csv", split_a[:radial_x])
NaCsData.dump_raw("$(sorted_prefix)_ry.csv", split_a[:radial_y])
NaCsData.dump_raw("$(sorted_prefix)_az.csv", split_a[:axial_z])

figure()
NaCsPlot.plot_survival_data(data_rx, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis X (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_ry, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Y (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(data_az, fmt="C0o-")
grid()
ylim([0, 0.5])
title("Axis Z (axial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()

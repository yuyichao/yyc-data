#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
import NaCsCalc.Format: Unc, Sci

const iname780 = joinpath(@__DIR__, "data", "aosense_780_20210706.csv")
const data780 = readdlm(iname780, ',', Float64, skipstart=1)
const iname935 = joinpath(@__DIR__, "data", "aosense_935_20210706.csv")
const data935 = readdlm(iname935, ',', Float64, skipstart=1)
const iname369 = joinpath(@__DIR__, "data", "aosense_369_20210706.csv")
const data369 = readdlm(iname369, ',', Float64, skipstart=1)
const iname650 = joinpath(@__DIR__, "data", "aosense_650_20210706.csv")
const data650 = readdlm(iname650, ',', Float64, skipstart=1)
const iname493 = joinpath(@__DIR__, "data", "aosense_493_20210706.csv")
const data493 = readdlm(iname493, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "aosense_20200806")

figure()
plot(data780[:, 1], data780[:, 2], "C0o-")
axvline(70, color="C3", ls="--")
text(70, 4.2, "Expected Threshold", rotation=90, va="center", ha="right",
     fontsize=11, color="C3")
grid()
ylim([0, 14.8])
xlabel("Current (mA)")
ylabel("Output Power (mW)")
title("Yb 780nm")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_780")

figure()
plot(data935[:, 1], data935[:, 2], "C0o-", label="Measure")
plot(88.6, 27.3, "C3X", label="Spec", markersize=15)
legend()
grid()
ylim([0, 29.5])
xlabel("Current (mA)")
ylabel("Output Power (mW)")
title("Yb 935nm")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_935")

figure()
plot(data369[:, 1], data369[:, 2], "C0o-", label="Measure")
plot(57.17, 2.94, "C3X", label="Spec", markersize=15)
legend()
grid()
ylim([0, 3.5])
xlabel("Current (mA)")
ylabel("Output Power (mW)")
title("Yb 369nm")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_369")

figure()
plot(data650[:, 1], data650[:, 2], "C0o-", label="Measure")
plot(130, 9.6, "C3X", label="Spec", markersize=15)
legend()
grid()
ylim([0, 13.9])
xlabel("Current (mA)")
ylabel("Output Power (mW)")
title("Ba 650nm")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_650")

figure()
plot(data493[:, 1], data493[:, 2], "C0o-", label="Measure")
plot(83.1, 12, "C3X", label="Spec", markersize=15)
legend()
grid()
ylim([0, 14.7])
xlabel("Current (mA)")
ylabel("Output Power (mW)")
title("Ba 493nm")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_493")

NaCsPlot.maybe_show()

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
using Printf
using SpecialFunctions

const name_profile1, (data_profile1,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000018906-RamanRabiCenterScan.h5"),
                            thresh=31)
const name_profile2, (data_profile2,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000018905-RamanRabiCenterScan.h5"),
                            thresh=31)
const name_profile3, (data_profile3,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000018934-RamanRabiCenterScan.h5"),
                            thresh=31)

const mid_pos = (data_profile1.params[1] + data_profile1.params[end]) / 2

const prefix = joinpath(@__DIR__, "imgs", "data_20230523_profile")

const data_profile1_shift = NaCsData.map_params((i, v)->v - mid_pos, data_profile1)
const data_profile2_shift = NaCsData.map_params((i, v)->v - mid_pos, data_profile2)
const data_profile3_shift = NaCsData.map_params((i, v)->v - mid_pos, data_profile3)

function line_profile_rabi(x, w, a0, lo)
    a = a0 * exp(-(x / w)^2)
    return lo + (1 - lo) * sin(a)^2
end

function model_line_profile(x, p)
    return line_profile_rabi.(x .- p[1], p[2], p[3], p[4])
end

const fit_profile1 = fit_loading(model_line_profile, data_profile1_shift,
                                 [-5, 1, π / 2, 0.05])
const fit_profile2 = fit_loading(model_line_profile, data_profile2_shift,
                                 [0, 1, π / 2, 0.05])
const fit_profile3 = fit_loading(model_line_profile, data_profile3_shift,
                                 [5, 1, π / 2, 0.05])
@show fit_profile1.uncs
@show fit_profile2.uncs
@show fit_profile3.uncs

figure()
NaCsPlot.plot_loading_data(data_profile1_shift, fmt="C0o")
plot(fit_profile1.plotx, fit_profile1.ploty, color="C0", label="beam 1")
NaCsPlot.plot_loading_data(data_profile2_shift, fmt="C1o")
plot(fit_profile2.plotx, fit_profile2.ploty, color="C1", label="beam 2")
NaCsPlot.plot_loading_data(data_profile3_shift, fmt="C2o")
plot(fit_profile3.plotx, fit_profile3.ploty, color="C2", label="beam 3")
grid()
legend(ncol=3, fontsize=12, loc="upper center")
ylim([0, 0.54])
xlim([-10.7, 10.7])
xlabel("Ion position (\$\\mu m\$)")
ylabel("\$P_{1}\$")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

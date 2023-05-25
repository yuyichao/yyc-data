#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Trap
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_loading
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using Printf
using SpecialFunctions

const name_435_rabi, (data_435_rabi,) =
    NaCsData.load_dax_scan1(joinpath(@__DIR__, "data/000006520-Rabi435.h5"))

const prefix = joinpath(@__DIR__, "imgs", "data_20221005_435_temp")

function get_model(ηs, Δns, fs)
    kernel(t, T) = Trap.thermal_sideband(max(T, 0.0) ./ fs, t, ηs, Δns)
    function model(x, p)
        return p[1] .+ p[2] .* kernel.((x .* 1e6 .- p[4]) ./ p[3], p[5])
    end
end

const N_A = 6.02214076e23
const m_Yb171 = 170.9363315e-3 / N_A # kg
const f_axial = 1.03e6
const f_radial1 = 1.52e6
const f_radial2 = 1.76e6
const f_all = (f_axial, f_radial1, f_radial2)
const λ_435 = 435e-9
const k_435 = 2π / λ_435

const Δns = (0, 0, 0)
const ηs = Trap.η.(m_Yb171, f_all, k_435 .* (sqrt(0.5), 0.5, 0.5))

const model_0 = get_model(ηs, (0, 0, 0), f_all ./ 1e6)
const model_ax = get_model(ηs[1], 0, f_all[1] ./ 1e6)

const fit_435_rabi = fit_loading(model_0, data_435_rabi,
                                 [1, -0.9, 5.5, 0, 20])
const fit_435_rabi_ax = fit_loading(model_ax, data_435_rabi,
                                    [1, -0.9, 5.5, 0, 20])

@show fit_435_rabi.uncs
@show fit_435_rabi_ax.uncs

figure()
NaCsPlot.plot_loading_data(data_435_rabi, xscale=1e6, fmt="C0o")
plot(fit_435_rabi.plotx .* 1e6, fit_435_rabi.ploty, color="C0", label="3D")
plot(fit_435_rabi_ax.plotx .* 1e6, fit_435_rabi_ax.ploty, color="C1", label="1D")
text(12.5, 0.02, "T\$_{3D}=$(fit_435_rabi.uncs[5]) MHz\$", color="C0")
text(12.5, 0.115, "T\$_{1D}=$(fit_435_rabi_ax.uncs[5]) MHz\$", color="C1")
grid()
ylim([0, 1])
xlabel("435 Pulse Time (\$\\mu s\$)")
ylabel("\$P_{detect}\$")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("../nacs_202003/molecular_raman_model.jl")

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200619_210941.mat",
                "data_20200621_001933.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [([0.04, 0.07, 0.085, 0.10, 0.13, 0.145, 0.16, 0.175, 0.19, 0.22, 0.25,
                 0.28] * 6 .- 0.01, # 3 mW, 770.3472 MHz
                [0.0], # 3 mW, 0 ms
                347.2 .+ (-2:0.4:2), # 3 mW, 1.02 ms
                347.2 .+ (-2:0.4:2), # 3 mW, 0.51 ms
                ),
               ([0.0, 2, 4, 10],
                [0.0, 40, 80, 120],
                [0.0, 40, 80, 160] .* 2,
                [0.0, 25, 50, 100],
                [0.0, 50, 100, 200],
                [0.0, 50, 100, 200] .* 2),
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((-1, 2), (-3, 4,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1, -2), (3, -4,)), maxcnts, specs)

const data_nacs_t = datas_nacs[1][1] # Survival 1
const data_nacs_00 = datas_nacs[1][2] # Survival 1
const data_nacs_101 = datas_nacs[1][3] # Survival 1
const data_nacs_50 = datas_nacs[1][4] # Survival 1

const data_pa_na = [datas_na[2][6]; datas_na[2][3]]
const data_pa_cs = datas_cs[2][6]
const data_pa_a = datas_nacs[2][6]
const data_pa_m = datas_nacs[2][3]

const data_fit = [NaCsData.map_params((i, v) -> (1, 347.2, v, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (1, v, 0.50, 1), data_nacs_50);
                  NaCsData.map_params((i, v) -> (1, v, 1.01, 1), data_nacs_101);
                  NaCsData.map_params((i, v) -> (1, 347.2, v, 1), data_nacs_t);
                  NaCsData.map_params((i, v) -> (2, v, 0.0, 2), data_pa_na);
                  NaCsData.map_params((i, v) -> (3, v, 0.0, 3), data_pa_cs);
                  NaCsData.map_params((i, v) -> (4, v, 0.0, 4), data_pa_a);
                  NaCsData.map_params((i, v) -> (5, v, 0.0, 5), data_pa_m)]

const prefix = joinpath(@__DIR__, "imgs", "fit_20200619_210941_raman_3322_2")

function model_exp(x, p)
    if length(p) > 2
        p[1] .* exp.(.- x .* p[2]) .+ p[3]
    else
        p[1] .* exp.(.- x .* p[2])
    end
end

function get_model_param(p, idx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[7 + idx]
    p0 = p1 * p0r
    return (p0, p1, f0, Ω, Γ1 + Γ_na + Γ_cs, Γ2)
end

function gen_pa_param(p, typ, pidx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[7 + pidx]
    if typ == 2
        return (p1, Γ_na)
    elseif typ == 3
        return (p1, Γ_cs)
    elseif typ == 4
        return (p1, Γ_na + Γ_cs)
    elseif typ == 5
        return (p1, Γ1 + Γ_na + Γ_cs)
    end
end

function model(xs, p)
    # p: f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r, p1s...
    function wrapper(x)
        typ = x[1]
        if typ == 1
            _, f, t, idx = x
            return model_2d(t, f, get_model_param(p, idx))
        else
            _, t, _, idx = x
            return model_exp(t, gen_pa_param(p, typ, idx))
        end
    end
    return wrapper.(xs)
end

fit = fit_survival(model, data_fit, [347.2, 2π * 1.0, 0, 2π * 0.3, 0, 0, 0.2,
                                     0.30, 0.8, 0.3, 0.3, 0.3],
                   plotx=false, lower=zeros(12))
@show fit.uncs
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)

const plot_freq_lo = 347.2 - 2.1
const plot_freq_hi = 347.2 + 2.1
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)

figure(figsize=[12.6, 11.2])

subplot(2, 2, 1)
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)
errorbar([param_1[3]], [ratio_00], [uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_50, fmt="C1.", label="0.50 ms")
plot(plot_freq, model_2d.(0.50, plot_freq, (param_1,)), "C1")
NaCsPlot.plot_survival_data(data_nacs_101, fmt="C2.", label="1.01 ms")
plot(plot_freq, model_2d.(1.01, plot_freq, (param_1,)), "C2")
legend(fontsize="x-small", loc="lower left")
xlim([plot_freq_lo, plot_freq_hi])
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288625 GHz, 3 mW")

subplot(2, 2, 2)
const plot_time = linspace(0, 1.70, 1000)
NaCsPlot.plot_survival_data([data_nacs_00; data_nacs_t], fmt="C0.", label="770.3472 MHz")
plot(plot_time, model_2d.(plot_time, 347.2, (param_1,)), "C0")
title("Time")
text(0.35, 0.20, ("\$f_{res}=$(uncs_1[3] / 1000 + 770)\$ MHz\n" *
                  "\$\\Omega_{Raman}=2\\pi\\times$(uncs_1[4] / 2π) kHz\$\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(uncs_1[5] / 2π * 1000)\$ Hz\n" *
                  "\$\\Gamma_{molecule}=2\\pi\\times$(uncs_1[6] / 2π) kHz\$\n"),
     color="C0", fontsize="small")
xlim([0, 1.75])
legend(fontsize="x-small", loc="upper right")
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.01\\ ms\\ offset)}}\$")
ylabel("Two-body survival")

subplot(2, 2, 3)
const plot_pa_time = linspace(0, 320, 1000)
NaCsPlot.plot_survival_data(data_pa_na, fmt="C0s", label="Na")
plot(plot_pa_time, model_exp(plot_pa_time, gen_pa_param(fit.param, 2, 2)), "C0-")
NaCsPlot.plot_survival_data(data_pa_cs, fmt="C1s", label="Cs")
plot(plot_pa_time, model_exp(plot_pa_time, gen_pa_param(fit.param, 3, 3)), "C1-")
NaCsPlot.plot_survival_data(data_pa_a, fmt="C2s", label="Hot")
plot(plot_pa_time, model_exp(plot_pa_time, gen_pa_param(fit.param, 4, 4)), "C2-")
NaCsPlot.plot_survival_data(data_pa_m, fmt="C3s", label="Cold")
plot(plot_pa_time, model_exp(plot_pa_time, gen_pa_param(fit.param, 5, 5)), "C3-")
xlim([0, 330])
ylim([0, 0.9])
legend(fontsize="small", ncol=2, borderpad=0.2, labelspacing=0.2,
       handletextpad=0.3, columnspacing=0.2, borderaxespad=0.4, loc="center left")
title("Atomic Lifetime")
grid()
xlabel("Time (ms)")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

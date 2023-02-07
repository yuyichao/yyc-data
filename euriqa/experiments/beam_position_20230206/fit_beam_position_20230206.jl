#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data
using NaCsPlot
using PyPlot
import NaCsCalc.Format: Unc, Sci

function load_positions(name)
    return readdlm(joinpath(@__DIR__, "data", "$(name).csv"),
                   ',', Float64, skipstart=1)
end

const pos_m1 = load_positions("m1_20230206")
const pos_m3 = load_positions("m3_20230206_2")

const npos_m1 = size(pos_m1, 1)
const npos_m3 = size(pos_m3, 1)

const center_x0 = 1430.0

function model_single(x, v, y0, y0′, b, b′)
    y0 = y0 + y0′ * v
    b = b + b′ * v
    return (x - center_x0) * b + y0
end

# Parameters:
# M1: α_pzt, y0, y0′, b, b′
# M3: y0, y0′, b, b′

function model_single(idx, α_pzt, y0_1, y0′_1, b_1, b′_1, y0_3, y0′_3, b_3, b′_3)
    if idx <= npos_m1
        return model_single(pos_m1[idx, 3], pos_m1[idx, 1] * α_pzt + pos_m1[idx, 2],
                            y0_1, y0′_1, b_1, b′_1) - pos_m1[idx, 4]
    end
    idx -= npos_m1
    if idx <= npos_m3
        return model_single(pos_m3[idx, 2], pos_m3[idx, 1],
                            y0_3, y0′_3, b_3, b′_3) - pos_m3[idx, 3]
    end
    idx -= npos_m3
    error("Index out of bound")
end

function model(idx, p)
    α_pzt, y0_1, y0′_1, b_1, b′_1, y0_3, y0′_3, b_3, b′_3 = p
    return model_single.(idx, α_pzt, y0_1, y0′_1, b_1, b′_1, y0_3, y0′_3, b_3, b′_3)
end

fit = fit_data(model, 1:(npos_m1 + npos_m3), zeros(npos_m1 + npos_m3),
               zeros(9), plotx=false)

function div_unc(a, b, cov)
    c = a / b
    v = [1 / b, -a / b^2]
    return Unc(c, sqrt(v' * cov * v))
end

println("""
α_pzt = $(fit.uncs[1])
y0_1′ = $(fit.uncs[3])
b_1′ = $(fit.uncs[5])
b_1′ / y0_1′ = $(div_unc(fit.param[5], fit.param[3], fit.covar[[5, 3], [5, 3]]))

y0_3′ = $(fit.uncs[7])
b_3′ = $(fit.uncs[9])
b_3′ / y0_3′ = $(div_unc(fit.param[9], fit.param[7], fit.covar[[9, 7], [9, 7]]))
""")

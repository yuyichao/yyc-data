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

const prefix = joinpath(@__DIR__, "imgs/ind_position_20230206")

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

function model_single_point(idx, α_pzt, y0_1, y0′_1, b_1, b′_1, y0_3, y0′_3, b_3, b′_3)
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
    return model_single_point.(idx, α_pzt, y0_1, y0′_1, b_1, b′_1,
                               y0_3, y0′_3, b_3, b′_3)
end

fit = fit_data(model, 1:(npos_m1 + npos_m3), zeros(npos_m1 + npos_m3),
               zeros(9), plotx=false)
fit_params = (α_pzt=fit.param[1],
              y0_1=fit.param[2], y0′_1=fit.param[3],
              b_1=fit.param[4], b′_1=fit.param[5],
              y0_3=fit.param[6], y0′_3=fit.param[7],
              b_3=fit.param[8], b′_3=fit.param[9])

function div_unc(a, b, cov)
    c = a / b
    v = [1 / b, -a / b^2]
    return Unc(c, sqrt(v' * cov * v))
end

print("""
α_pzt = $(fit.uncs[1])
y0_1′ = $(fit.uncs[3])
b_1′ = $(fit.uncs[5])
b_1′ / y0_1′ = $(div_unc(fit.param[5], fit.param[3], fit.covar[[5, 3], [5, 3]]))

y0_3′ = $(fit.uncs[7])
b_3′ = $(fit.uncs[9])
b_3′ / y0_3′ = $(div_unc(fit.param[9], fit.param[7], fit.covar[[9, 7], [9, 7]]))
""")

function group_data(pos)
    data = Dict{Vector{Float64},NTuple{2,Vector{Float64}}}()
    for i in 1:size(pos, 1)
        key = pos[i, 1:end - 2]
        x = pos[i, end - 1]
        y = pos[i, end]
        value = get!(()->(Float64[], Float64[]), data, key)
        push!(value[1], x)
        push!(value[2], y)
    end
    return sort!([(k, v[1], v[2]) for (k, v) in data])
end

lines_m1 = group_data(pos_m1)
lines_m3 = group_data(pos_m3)

const plotx_1 = range(830, 2030, length=101)
figure()
for i in 1:length(lines_m1)
    k, x, y = lines_m1[i]
    plot(x, y, "C$(i)o", label="$(Int(k[1])) V, $(Int(k[2])) Turn")
    plot(plotx_1, model_single.(plotx_1, k[1] * fit_params.α_pzt + k[2],
                                fit_params.y0_1, fit_params.y0′_1,
                                fit_params.b_1, fit_params.b′_1), "C$(i)")
end
grid()
legend(ncol=2, fontsize=12)
NaCsPlot.maybe_save("$(prefix)_m1")

const plotx_3 = range(170, 1550, length=101)
figure()
for i in 1:length(lines_m3)
    k, x, y = lines_m3[i]
    plot(x, y, "C$(i)o", label="$(k[1]) deg")
    plot(plotx_3, model_single.(plotx_3, k[1], fit_params.y0_3, fit_params.y0′_3,
                                fit_params.b_3, fit_params.b′_3), "C$(i)")
end
grid()
legend(ncol=2, fontsize=12)
NaCsPlot.maybe_save("$(prefix)_m3")

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("read_data.jl")

using NaCsData
using NaCsPlot
using PyPlot

const prefix = joinpath(@__DIR__, "imgs/scope_trace_20230925")

const data1 = read_bin_compressed(joinpath(@__DIR__, "data", "cavity_20230925_1_1.zst"))
const data2 = read_bin_compressed(joinpath(@__DIR__, "data", "cavity_20230925_2_1.zst"))
const data3_1 = read_bin_compressed(joinpath(@__DIR__, "data",
                                             "cavity_20230925_3_1.zst"))
const data3_2 = read_bin_compressed(joinpath(@__DIR__, "data",
                                             "cavity_20230925_3_2.zst"))
const data4_1 = read_bin_compressed(joinpath(@__DIR__, "data",
                                             "cavity_20230925_4_1.zst"))
const data4_2 = read_bin_compressed(joinpath(@__DIR__, "data",
                                             "cavity_20230925_4_2.zst"))

function get_xy_data(data; xscale=1, yscale=1)
    xs = ((0:(length(data.data) - 1)) .* data.dx .+ data.x0) .* xscale
    ys = (data.data .* data.dy .+ data.y0) .* yscale
    return xs, ys
end

const xy1 = get_xy_data(data1, xscale=1000)
const xy2 = get_xy_data(data2, xscale=1000)
const xy3_1 = get_xy_data(data3_1, xscale=1000)
const xy3_2 = get_xy_data(data3_2, xscale=1000)
const xy4_1 = get_xy_data(data4_1, xscale=1000)
const xy4_2 = get_xy_data(data4_2, xscale=1000)

figure(figsize=[6.4 * 100, 4.8])
plot(xy1..., color="C0")
ylim([0, 0.065])
ylabel("Transmission")
grid()
xlim([minimum(xy1[1]), maximum(xy1[1])])
xlabel("Time (ms)")
NaCsPlot.maybe_save("$(prefix)_1")

figure(figsize=[6.4 * 100, 4.8])
plot(xy2..., color="C0")
ylim([0, 0.076])
ylabel("Transmission")
grid()
xlim([minimum(xy2[1]), maximum(xy2[1])])
xlabel("Time (ms)")
NaCsPlot.maybe_save("$(prefix)_2")

figure(figsize=[6.4 * 100, 4.8])
plot(xy3_2..., color="C0")
ylim([0.42, 0.565])
ylabel("Reflection", color="C0")
grid()
ax2 = gca()[:twinx]()
plot(xy3_1..., color="C1")
ylim([0, 0.0725])
ylabel("Transmission", color="C1")
xlim([minimum(xy3_1[1]), maximum(xy3_1[1])])
xlabel("Time (ms)")
NaCsPlot.maybe_save("$(prefix)_3")

figure(figsize=[6.4 * 3, 4.8])
plot(xy4_2..., color="C0")
ylim([0.44, 0.57])
ylabel("Reflection", color="C0")
grid()
ax2 = gca()[:twinx]()
plot(xy4_1..., color="C1")
ylim([0, 0.065])
ylabel("Transmission", color="C1")
xlim([-0.2, 0.3])
xlabel("Time (ms)")
NaCsPlot.maybe_save("$(prefix)_4")

NaCsPlot.maybe_show()

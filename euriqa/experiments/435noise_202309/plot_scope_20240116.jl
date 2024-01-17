#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("read_data.jl")

using NaCsData
using NaCsPlot
using PyPlot

const prefix = joinpath(@__DIR__, "imgs/scope_trace_20240116")

const data1 = read_bin_compressed(joinpath(@__DIR__, "data", "cavity_20240116_5Hz_signal_only_1.zst"))
const data2 = read_bin_compressed(joinpath(@__DIR__, "data", "cavity_20240116_5Hz_signal_only_2.zst"))

function get_xy_data(data; xscale=1, yscale=1)
    xs = ((0:(length(data.data) - 1)) .* data.dx .+ data.x0) .* xscale
    ys = (data.data .* data.dy .+ data.y0) .* yscale
    return xs, ys
end

const xy1 = get_xy_data(data1, xscale=1000)
const xy2 = get_xy_data(data2, xscale=1000)

figure(figsize=[6.4 * 100, 4.8])
plot(xy1..., color="C0")
plot(xy2..., color="C1")
ylabel("Transmission")
grid()
xlim([minimum(xy1[1]), maximum(xy1[1])])
xlabel("Time (ms)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

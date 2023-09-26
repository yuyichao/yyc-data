#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("read_data.jl")

using NaCsData
using NaCsPlot
using PyPlot
using Statistics

const prefix = joinpath(@__DIR__, "imgs/scope_trace_822_20230926")

const data1 = read_bin_compressed(joinpath(@__DIR__, "data",
                                           "822_cavity_20230926_1_1.zst"))
const data2 = read_bin_compressed(joinpath(@__DIR__, "data",
                                           "822_cavity_20230926_1_2.zst"))

function get_xy_data(data; xscale=1, yscale=1)
    xs = ((0:(length(data.data) - 1)) .* data.dx .+ data.x0) .* xscale
    ys = (data.data .* data.dy .+ data.y0) .* yscale
    return xs, ys
end

moving_average(data, window) =
    [mean(@view(data[i:i + window - 1])) for i in 1:window:(length(data) - window + 1)]

const xy1 = get_xy_data(data1, xscale=1000)
const xy2 = get_xy_data(data2, xscale=1000)

figure()
plot(moving_average(xy1[1], 100),
     moving_average(xy1[2], 100), color="C0")
ylim([-0.08, 0.08])
ylabel("Error Signal", color="C0")
grid()
xlabel("Time (ms)")
ax2 = gca()[:twinx]()
plot(xy2..., color="C1")
ylim([0, 2.0])
ylabel("Transmission", color="C1")
xlim([minimum(xy1[1]), maximum(xy1[1])])
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot
using DataStructures

matplotlib["rcParams"][:update](Dict("font.size" => 20))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

const dts = [36, 132, 84, 64, 69, 87, 2, 13, 70, 17, 68, 10, 41, 36, 2,
             98, 17, 33, 190, 14, 69, 32, 169, 46, 28, 15, 34, 41]

function get_plot_data(dts)
    t = 0
    p = 1
    ts = Float64[t]
    ps = Int[p]
    for dt in dts
        t += dt
        push!(ts, t)
        push!(ts, t)
        push!(ps, p)
        push!(ps, -p)
        p = -p
    end
    return ts, ps
end

const ts, ps = get_plot_data(dts)

const prefix = joinpath(@__DIR__, "../imgs/double-well")

figure()
plot(ts, ps)
grid(axis="x")
yticks([-1, 1])
ylim([-1.2, 1.2])
xlabel("Time (min)")
ylabel("Ion Position")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

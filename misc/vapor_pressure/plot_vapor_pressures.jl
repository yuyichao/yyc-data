#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("vapor_data.jl")

using PyPlot
using NaCsPlot

const prefix = joinpath(@__DIR__, "img/vapor_pressures")

function plot_data(Tlo, Thi, data, color, label)
    if data.Tmelt > Thi
        Ts = range(Tlo, Thi, 10000)
        plot(Ts, vapor_pressure.(Ts, Ref(data)), color=color, label=label)
    else
        Ts = range(Tlo, data.Tmelt, 10000)
        plot(Ts, vapor_pressure.(Ts, Ref(data)), color=color, label=label)
        Ts = range(data.Tmelt, Thi, 10000)
        plot(Ts, vapor_pressure.(Ts, Ref(data)), "--", color=color)
        plot(data.Tmelt, vapor_pressure(data.Tmelt, data), "o", color=color)
    end
end

figure()
gca().add_patch(matplotlib.patches.Rectangle((50, 5e-5), 860 - 50, 2e-3 - 5e-5,
                                             fill=true, alpha=0.2, facecolor="g",
                                             edgecolor=nothing))
plot_data(60, 850, NaData, "C0", "Na")
plot_data(60, 850, YbData, "C1", "Yb")
plot_data(60, 850, LiData, "C2", "Li")
plot_data(60, 850, CaData, "C3", "Ca")
plot_data(60, 850, TlData, "C4", "Tl")
plot_data(60, 850, AgData, "C5", "Ag")
plot_data(60, 850, InData, "C6", "In")
legend(fontsize=12)
yscale("log")
grid()
xlabel("Temperature (\$\\!\\!^\\circ\\!\\!\$C)")
ylabel("Vapor pressure (Torr)")
xlim([60, 850])
ylim([5e-12, 1.5])
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

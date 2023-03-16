#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc: Trap
using NaCsPlot
using PyPlot

matplotlib["rcParams"][:update](Dict("font.weight" => "normal"))

const m_Yb = 171e-3 / 6.02e23
const ω_m = 2π * 200e3
const η = Trap.η(m_Yb, ω_m / (2π), 2π / 370e-9)

function doppler_spectrum(Δn, nbar, η)
    return exp(-(Δn - η^2)^2 / 4 / η^2 / nbar) / 2 / η / sqrt(π * nbar)
end

function sideband_spectrum(Δn, nbar, η)
    nmax = round(Int, nbar * 20 + 5)
    r = 0.0
    for n in 0:nmax
        r += Trap.sideband(n, n + Δn, η)^2 * exp(-n / nbar)
    end
    return r * (1 - exp(-1 / nbar))
end

const Δns = -20:20
const Δns_cont = range(-20, 20, length=2001)

figure()
plot(Δns_cont, doppler_spectrum.(Δns_cont, 0.5, η), "C0")
plot(Δns, sideband_spectrum.(Δns, 0.5, η), "C0o", label="\$\\bar n=0.5\$")

plot(Δns_cont, doppler_spectrum.(Δns_cont, 5, η) .* sqrt(5), "C1")
plot(Δns, sideband_spectrum.(Δns, 5, η) .* sqrt(5), "C1o", label="\$\\bar n=5\$")

plot(Δns_cont, doppler_spectrum.(Δns_cont, 50, η) .* sqrt(50), "C2")
plot(Δns, sideband_spectrum.(Δns, 50, η) .* sqrt(50), "C2o", label="\$\\bar n=50\$")

plot(Δns_cont, doppler_spectrum.(Δns_cont, 500, η) .* sqrt(500), "C3")
plot(Δns, sideband_spectrum.(Δns, 500, η) .* sqrt(500), "C3o", label="\$\\bar n=500\$")
legend(fontsize=13, ncol=2)
grid()
xlabel("Detuning")
ylabel("Coupling strength")
NaCsPlot.maybe_save(joinpath(@__DIR__, "../imgs/sideband_vs_doppler"))

function plot_sidebands(ns, Δns, η)
    for Δn in Δns
        plot(ns, abs.(Trap.sideband.(ns, ns .+ Δn, η)), ".-",
             label="\$\\Delta n=$(Δn)\$")
    end
end
figure()
plot_sidebands(0:100, 0:-1:-3, η)
xlim([0, 100])
ylim([0, 1])
grid()
ylabel("\$|\\langle n |e^{ik\\hat x}| n + \\Delta n \\rangle|\$")
text(3, 0.88, "\$\\Delta n\\!\\!=\\!0\$", color="C0")
text(12, 0.61, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!1\$", color="C1")
text(47, 0.51, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!2\$", color="C2")
text(80, 0.46, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!3\$", color="C3")
gca()[:get_yaxis]()[:set_label_coords](-0.105, 0.5)
xlabel("Motional state \$n\$")
NaCsPlot.maybe_save(joinpath(@__DIR__, "../imgs/sideband_coupling"))

NaCsPlot.maybe_show()

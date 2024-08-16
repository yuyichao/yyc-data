#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

# Unit: us, MHz
function get_sys(Ω, Δ, B)
    sys = System{ComplexF64}(4)

    H = new_H(sys)
    # S1/2: (1, -1), (1, 0), (1, 1)
    # Detuning
    H_add_E!(H, Δ, 1)
    H_add_E!(H, Δ, 2)
    H_add_E!(H, Δ, 3)
    # Zeeman
    H_add_E!(H, 2π * -1.4 * B, 1)
    H_add_E!(H, 2π * 1.4 * B, 3)

    Ωπ = Ω / sqrt(2)
    Ωσ = Ω / 2
    # π transition
    # (1, 0) -> (0, 0)
    H_add_Ω!(H, Ωπ, 2, 4)

    # σ⁺ transition
    # (1, -1) -> (0, 0)
    H_add_Ω!(H, Ωσ, 1, 4)

    # σ⁻ transition
    # (1, 1) -> (0, 0)
    H_add_Ω!(H, Ωσ, 3, 4)
    add_H!(sys, H)

    Γ = 2π * 19.6
    # π transition
    Cπ = new_C(sys)
    # (1, 0) -> (0, 0)
    C_add_decay!(Cπ, Γ / 3, 2, 4)
    add_C!(sys, Cπ)

    # σ⁺ transition
    Cσ⁺ = new_C(sys)
    # (1, -1) -> (0, 0)
    C_add_decay!(Cσ⁺, Γ / 3, 1, 4)
    add_C!(sys, Cσ⁺)

    # σ⁻ transition
    Cσ⁻ = new_C(sys)
    # (1, 1) -> (0, 0)
    C_add_decay!(Cσ⁻, Γ / 3, 3, 4)
    add_C!(sys, Cσ⁻)

    return sys
end

const ρ0 = zeros(4, 4)
ρ0[2, 2] = 1

function get_pe(Ω, B)
    sys = get_sys(Ω, 0, B)
    ρ = propagate(sys, ρ0, 3.0)
    return real(ρ[4, 4])
end

const Bs = range(0, 20, 101)
const Ωs = range(0, 2π * 60, 101)

pes = [get_pe(Ω, B) for Ω in Ωs, B in Bs]

figure()
imshow(pes, origin="lower", aspect="auto",
       extent=(-step(Bs) / 2, last(Bs) + step(Bs) / 2,
               -step(Ωs) / 2 / 2π, (last(Ωs) + step(Ωs) / 2) / 2π))
xlabel("B (G)")
ylabel("\$\\Omega/2\\pi (MHz)\$")
colorbar()

# const sys_1 = get_sys(2π * 1, 0, B)
# const sys_2 = get_sys(2π * 2, 0, B)
# const sys_5 = get_sys(2π * 5, 0, B)
# const sys_10 = get_sys(2π * 10, 0, B)
# const sys_20 = get_sys(2π * 20, 0, B)
# const sys_50 = get_sys(2π * 50, 0, B)
# const sys_100 = get_sys(2π * 100, 0, B)
# const sys_200 = get_sys(2π * 200, 0, B)
# const sys_500 = get_sys(2π * 500, 0, B)

# const ts = range(0, 2, 2001)

# function get_plot_data(sys, ρ0, ts)
#     p0 = Float64[]
#     p1 = Float64[]
#     pe = Float64[]
#     for t in ts
#         ρ = propagate(sys, ρ0, t)
#         push!(p0, 0)
#         push!(p1, real(ρ[1, 1]) + real(ρ[2, 2]) + real(ρ[3, 3]))
#         push!(pe, real(ρ[4, 4]))
#     end
#     return p0, p1, pe
# end

# const data_1 = @time get_plot_data(sys_1, ρ0, ts)
# const data_2 = @time get_plot_data(sys_2, ρ0, ts)
# const data_5 = @time get_plot_data(sys_5, ρ0, ts)
# const data_10 = @time get_plot_data(sys_10, ρ0, ts)
# const data_20 = @time get_plot_data(sys_20, ρ0, ts)
# const data_50 = @time get_plot_data(sys_50, ρ0, ts)
# const data_100 = @time get_plot_data(sys_100, ρ0, ts)
# const data_200 = @time get_plot_data(sys_200, ρ0, ts)
# const data_500 = @time get_plot_data(sys_500, ρ0, ts)

# # const prefix = joinpath(@__DIR__, "../imgs/three-states-scatter")

# # figure()
# # plot(ts, data_1[1], label="1")
# # plot(ts, data_2[1], label="2")
# # plot(ts, data_5[1], label="5")
# # plot(ts, data_10[1], label="10")
# # plot(ts, data_20[1], label="20")
# # plot(ts, data_50[1], label="50")
# # plot(ts, data_100[1], label="100")
# # plot(ts, data_200[1], label="200")
# # plot(ts, data_500[1], label="500")
# # ylim([0, 1])
# # grid()
# # xlabel("t")
# # ylabel("p0")
# # legend(fontsize=8)

# # figure()
# # plot(ts, data_1[2], label="1")
# # plot(ts, data_2[2], label="2")
# # plot(ts, data_5[2], label="5")
# # plot(ts, data_10[2], label="10")
# # plot(ts, data_20[2], label="20")
# # plot(ts, data_50[2], label="50")
# # plot(ts, data_100[2], label="100")
# # plot(ts, data_200[2], label="200")
# # plot(ts, data_500[2], label="500")
# # ylim([0, 1])
# # grid()
# # xlabel("t")
# # ylabel("p1")
# # legend(fontsize=8)

# figure()
# plot(ts, data_1[3], label="1")
# plot(ts, data_2[3], label="2")
# plot(ts, data_5[3], label="5")
# plot(ts, data_10[3], label="10")
# plot(ts, data_20[3], label="20")
# plot(ts, data_50[3], label="50")
# plot(ts, data_100[3], label="100")
# plot(ts, data_200[3], label="200")
# plot(ts, data_500[3], label="500")
# ylim([0, 1])
# grid()
# xlabel("t")
# ylabel("pe")
# legend(fontsize=8)
# # NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

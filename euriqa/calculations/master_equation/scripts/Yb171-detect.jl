#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

# Unit: us, MHz
function get_sys(Ω, Δ, B)
    sys = System{ComplexF64}(8)

    H = new_H(sys)
    # S1/2: (0, 0), (1, -1), (1, 0), (1, 1)
    # Detuning
    H_add_E!(H, Δ, 1)
    H_add_E!(H, Δ, 2)
    H_add_E!(H, Δ, 3)
    H_add_E!(H, Δ, 4)
    # HF
    H_add_E!(H, 2π * -12.643e3, 1)
    # Zeeman
    H_add_E!(H, 2π * -1.4 * B, 2)
    H_add_E!(H, 2π * 1.4 * B, 4)

    # P1/2: (0, 0), (1, -1), (1, 0), (1, 1)
    # HF
    H_add_E!(H, 2π * 2.105e3, 6)
    H_add_E!(H, 2π * 2.105e3, 7)
    H_add_E!(H, 2π * 2.105e3, 8)
    # Zeeman
    H_add_E!(H, 2π * -1.4 * B / 3, 6)
    H_add_E!(H, 2π * 1.4 * B / 3, 8)

    Ωπ = Ω / 2
    Ωσ = Ω / sqrt(2)
    # π transition
    # (0, 0) -> (1, 0)
    H_add_Ω!(H, Ωπ, 1, 7)
    # (1, -1) -> (1, -1)
    H_add_Ω!(H, -Ωπ, 2, 6)
    # (1, 0) -> (0, 0)
    H_add_Ω!(H, Ωπ, 3, 5)
    # (1, 1) -> (1, 1)
    H_add_Ω!(H, Ωπ, 4, 8)

    # σ⁺ transition
    # (0, 0) -> (1, 1)
    H_add_Ω!(H, Ωσ, 1, 8)
    # (1, -1) -> (0, 0)
    H_add_Ω!(H, Ωσ, 2, 5)
    # (1, -1) -> (1, 0)
    H_add_Ω!(H, Ωσ, 2, 7)
    # (1, 0) -> (1, 1)
    H_add_Ω!(H, -Ωσ, 3, 8)

    # σ⁻ transition
    # (0, 0) -> (1, -1)
    H_add_Ω!(H, Ωσ, 1, 6)
    # (1, 0) -> (1, -1)
    H_add_Ω!(H, Ωσ, 3, 6)
    # (1, 1) -> (0, 0)
    H_add_Ω!(H, Ωσ, 4, 5)
    # (1, 1) -> (1, 0)
    H_add_Ω!(H, -Ωσ, 4, 7)
    add_H!(sys, H)

    Γ = 2π * 19.6
    # π transition
    Cπ = new_C(sys)
    # (0, 0) -> (1, 0)
    C_add_decay!(Cπ, Γ / 3, 1, 7)
    # (1, -1) -> (1, -1)
    C_add_decay!(Cπ, Γ / 3, 2, 6)
    # (1, 0) -> (0, 0)
    C_add_decay!(Cπ, Γ / 3, 3, 5)
    # (1, 1) -> (1, 1)
    C_add_decay!(Cπ, Γ / 3, 4, 8)
    add_C!(sys, Cπ)

    # σ⁺ transition
    Cσ⁺ = new_C(sys)
    # (0, 0) -> (1, 1)
    C_add_decay!(Cσ⁺, Γ / 3, 1, 8)
    # (1, -1) -> (0, 0)
    C_add_decay!(Cσ⁺, Γ / 3, 2, 5)
    # (1, -1) -> (1, 0)
    C_add_decay!(Cσ⁺, Γ / 3, 2, 7)
    # (1, 0) -> (1, 1)
    C_add_decay!(Cσ⁺, Γ / 3, 3, 8)
    add_C!(sys, Cσ⁺)

    # σ⁻ transition
    Cσ⁻ = new_C(sys)
    # (0, 0) -> (1, -1)
    C_add_decay!(Cσ⁻, Γ / 3, 1, 6)
    # (1, 0) -> (1, -1)
    C_add_decay!(Cσ⁻, Γ / 3, 3, 6)
    # (1, 1) -> (0, 0)
    C_add_decay!(Cσ⁻, Γ / 3, 4, 5)
    # (1, 1) -> (1, 0)
    C_add_decay!(Cσ⁻, Γ / 3, 4, 7)
    add_C!(sys, Cσ⁻)

    return sys
end

const sys_1 = get_sys(2π * 1, 0, 3)
const sys_2 = get_sys(2π * 2, 0, 3)
const sys_5 = get_sys(2π * 5, 0, 3)
const sys_10 = get_sys(2π * 10, 0, 3)
const sys_20 = get_sys(2π * 20, 0, 3)
const sys_50 = get_sys(2π * 50, 0, 3)
const sys_100 = get_sys(2π * 100, 0, 3)
const sys_200 = get_sys(2π * 200, 0, 3)
const sys_500 = get_sys(2π * 500, 0, 3)

const ts = range(0, 2, 2001)
const ρ0 = zeros(8, 8)
ρ0[3, 3] = 1

function get_plot_data(sys, ρ0, ts)
    p0 = Float64[]
    p1 = Float64[]
    pe = Float64[]
    for t in ts
        ρ = propagate(sys, ρ0, t)
        push!(p0, real(ρ[1, 1]))
        push!(p1, real(ρ[2, 2]) + real(ρ[3, 3]) + real(ρ[4, 4]))
        push!(pe, real(ρ[5, 5]) + real(ρ[6, 6]) + real(ρ[7, 7]) + real(ρ[8, 8]))
    end
    return p0, p1, pe
end

const data_1 = @time get_plot_data(sys_1, ρ0, ts)
const data_2 = @time get_plot_data(sys_2, ρ0, ts)
const data_5 = @time get_plot_data(sys_5, ρ0, ts)
const data_10 = @time get_plot_data(sys_10, ρ0, ts)
const data_20 = @time get_plot_data(sys_20, ρ0, ts)
const data_50 = @time get_plot_data(sys_50, ρ0, ts)
const data_100 = @time get_plot_data(sys_100, ρ0, ts)
const data_200 = @time get_plot_data(sys_200, ρ0, ts)
const data_500 = @time get_plot_data(sys_500, ρ0, ts)

# const prefix = joinpath(@__DIR__, "../imgs/three-states-scatter")

figure()
plot(ts, data_1[1], label="1")
plot(ts, data_2[1], label="2")
plot(ts, data_5[1], label="5")
plot(ts, data_10[1], label="10")
plot(ts, data_20[1], label="20")
plot(ts, data_50[1], label="50")
plot(ts, data_100[1], label="100")
plot(ts, data_200[1], label="200")
plot(ts, data_500[1], label="500")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p0")
legend(fontsize=8)

figure()
plot(ts, data_1[2], label="1")
plot(ts, data_2[2], label="2")
plot(ts, data_5[2], label="5")
plot(ts, data_10[2], label="10")
plot(ts, data_20[2], label="20")
plot(ts, data_50[2], label="50")
plot(ts, data_100[2], label="100")
plot(ts, data_200[2], label="200")
plot(ts, data_500[2], label="500")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p1")
legend(fontsize=8)

figure()
plot(ts, data_1[3], label="1")
plot(ts, data_2[3], label="2")
plot(ts, data_5[3], label="5")
plot(ts, data_10[3], label="10")
plot(ts, data_20[3], label="20")
plot(ts, data_50[3], label="50")
plot(ts, data_100[3], label="100")
plot(ts, data_200[3], label="200")
plot(ts, data_500[3], label="500")
ylim([0, 1])
grid()
xlabel("t")
ylabel("pe")
legend(fontsize=8)
# NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

const HF_3D1 = 2π * 882
const g_3P0 = 2π * 1e-3 # ?
const g_3D1_1_2 = 2 / 3 * 2π * 1.4
const g_3D1_3_2 = 1 / 3 * 2π * 1.4

function E_3P0(B, mF2)
    return g_3P0 * B * mF2 / 2
end

function E_3D1(B, F2, mF2)
    if F2 == 1
        return g_3D1_1_2 * B * mF2 / 2
    else
        @assert F2 == 3
        return g_3D1_3_2 * B * mF2 / 2 + HF_3D1
    end
end

# Unit: us, MHz
function get_sys(Ω, Δ, B)
    # States:
    # 3P0 F1/2: -1/2, 1/2
    # 3D1 F1/2: -1/2, 1/2
    # 3D1 F3/2: -3/2, -1/2, 1/2, 3/2
    # Dummy state

    sys = System{ComplexF64}(9)

    Γe = 2π * 0.48 # 2 MHz
    branch_back = 0.64

    H = new_H(sys)
    # 3P0 F1/2: -1/2, 1/2
    H_add_E!(H, E_3P0(B, -1), 1)
    H_add_E!(H, E_3P0(B, 1), 2)
    # 3D1 F1/2: -1/2, 1/2
    H_add_E!(H, E_3D1(B, 1, -1) - Δ, 3)
    H_add_E!(H, E_3D1(B, 1, 1) - Δ, 4)
    # 3D1 F3/2: -3/2, -1/2, 1/2, 3/2
    H_add_E!(H, E_3D1(B, 3, -3) - Δ, 5)
    H_add_E!(H, E_3D1(B, 3, -1) - Δ, 6)
    H_add_E!(H, E_3D1(B, 3, 1) - Δ, 7)
    H_add_E!(H, E_3D1(B, 3, 3) - Δ, 8)
    # Dummy state

    Ωπ = Ω / sqrt(2) # * im
    Ωσ = Ω / 2
    # π transition
    # (1/2, -1/2) -> (1/2, -1/2)
    H_add_Ω!(H, sqrt(1 / 3) * Ωπ, 1, 3)
    # (1/2, 1/2) -> (1/2, 1/2)
    H_add_Ω!(H, -sqrt(1 / 3) * Ωπ, 2, 4)
    # (1/2, -1/2) -> (3/2, -1/2)
    H_add_Ω!(H, sqrt(2 / 3) * Ωπ, 1, 6)
    # (1/2, 1/2) -> (3/2, 1/2)
    H_add_Ω!(H, sqrt(2 / 3) * Ωπ, 2, 7)

    # σ⁺ transition
    # (1/2, -1/2) -> (1/2, 1/2)
    H_add_Ω!(H, sqrt(2 / 3) * Ωσ, 1, 4)
    # (1/2, -1/2) -> (3/2, 1/2)
    H_add_Ω!(H, sqrt(1 / 3) * Ωσ, 1, 7)
    # (1/2, 1/2) -> (3/2, 3/2)
    H_add_Ω!(H, Ωσ, 2, 8)

    # σ⁻ transition
    # (1/2, 1/2) -> (1/2, -1/2)
    H_add_Ω!(H, -sqrt(2 / 3) * Ωσ, 2, 3)
    # (1/2, -1/2) -> (3/2, -3/2)
    H_add_Ω!(H, Ωσ, 1, 5)
    # (1/2, 1/2) -> (3/2, -1/2)
    H_add_Ω!(H, sqrt(1 / 3) * Ωσ, 2, 6)
    add_H!(sys, H)

    # π transition
    Cπ = new_C(sys)
    # (1/2, -1/2) -> (1/2, -1/2)
    C_add_decay!(Cπ, (1 / 3) * Γe * branch_back, 1, 3)
    # (1/2, 1/2) -> (1/2, 1/2)
    C_add_decay!(Cπ, (1 / 3) * Γe * branch_back, 2, 4)
    # (1/2, -1/2) -> (3/2, -1/2)
    C_add_decay!(Cπ, (2 / 3) * Γe * branch_back, 1, 6)
    # (1/2, 1/2) -> (3/2, 1/2)
    C_add_decay!(Cπ, (2 / 3) * Γe * branch_back, 2, 7)
    add_C!(sys, Cπ)

    # σ⁺ transition
    Cσ⁺ = new_C(sys)
    # (1/2, -1/2) -> (1/2, 1/2)
    C_add_decay!(Cσ⁺, (2 / 3) * Γe * branch_back, 1, 4)
    # (1/2, -1/2) -> (3/2, 1/2)
    C_add_decay!(Cσ⁺, (1 / 3) * Γe * branch_back, 1, 7)
    # (1/2, 1/2) -> (3/2, 3/2)
    C_add_decay!(Cσ⁺, Γe * branch_back, 2, 8)
    add_C!(sys, Cσ⁺)

    # σ⁻ transition
    Cσ⁻ = new_C(sys)
    # (1/2, 1/2) -> (1/2, -1/2)
    C_add_decay!(Cσ⁻, (2 / 3) * Γe * branch_back, 2, 3)
    # (1/2, -1/2) -> (3/2, -3/2)
    C_add_decay!(Cσ⁻, Γe * branch_back, 1, 5)
    # (1/2, 1/2) -> (3/2, -1/2)
    C_add_decay!(Cσ⁻, (1 / 3) * Γe * branch_back, 2, 6)
    add_C!(sys, Cσ⁻)

    for i in 3:8
        C = new_C(sys)
        C_add_decay!(C, Γe * (1 - branch_back), 9, i)
        add_C!(sys, C)
    end

    return sys
end

const Ω_total = 2π * 2 / sqrt(2) * sqrt(100)
const B0 = 120

const sys_1_1 = get_sys(Ω_total, E_3D1(B0, 1, -1) - E_3P0(B0, -1), B0)
const sys_1_2 = get_sys(Ω_total, E_3D1(B0, 1, 1) - E_3P0(B0, -1), B0)

const sys_2_1 = get_sys(Ω_total, E_3D1(B0, 3, -3) - E_3P0(B0, -1), B0)
const sys_2_2 = get_sys(Ω_total, E_3D1(B0, 3, -1) - E_3P0(B0, -1), B0)
const sys_2_3 = get_sys(Ω_total, E_3D1(B0, 3, 1) - E_3P0(B0, -1), B0)
const sys_2_4 = get_sys(Ω_total, E_3D1(B0, 3, 3) - E_3P0(B0, -1), B0)

const ts = range(0, 1500, 1001)
const ρ0 = zeros(9, 9)
ρ0[2, 2] = 1

function get_plot_data(sys, ρ0, ts)
    p0 = Float64[]
    p1 = Float64[]
    pg = Float64[]
    for t in ts
        ρ = propagate(sys, ρ0, t)
        push!(p0, real(ρ[1, 1]))
        push!(p1, real(ρ[2, 2]))
        push!(pg, real(ρ[9, 9]) * 0.99) # account for leakage into 3P2
    end
    return p0, p1, pg
end

const data_1_1 = @time get_plot_data(sys_1_1, ρ0, ts)
const data_1_2 = @time get_plot_data(sys_1_2, ρ0, ts)

const data_2_1 = @time get_plot_data(sys_2_1, ρ0, ts)
const data_2_2 = @time get_plot_data(sys_2_2, ρ0, ts)
const data_2_3 = @time get_plot_data(sys_2_3, ρ0, ts)
const data_2_4 = @time get_plot_data(sys_2_4, ρ0, ts)

# const prefix = joinpath(@__DIR__, "../imgs/three-states-scatter")

figure()
plot(ts, data_1_1[1], label="1/2, -1/2")
plot(ts, data_1_2[1], label="1/2, +1/2")
plot(ts, data_2_1[1], label="3/2, -3/2")
plot(ts, data_2_2[1], label="3/2, -1/2")
plot(ts, data_2_3[1], label="3/2, +1/2")
plot(ts, data_2_4[1], label="3/2, +3/2")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p(P,-1/2)")
legend(fontsize=8)

figure()
plot(ts, data_1_1[2], label="1/2, -1/2")
plot(ts, data_1_2[2], label="1/2, +1/2")
plot(ts, data_2_1[2], label="3/2, -3/2")
plot(ts, data_2_2[2], label="3/2, -1/2")
plot(ts, data_2_3[2], label="3/2, +1/2")
plot(ts, data_2_4[2], label="3/2, +3/2")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p(P,+1/2)")
legend(fontsize=8)

figure()
plot(ts, data_1_1[3], label="1/2, -1/2")
plot(ts, data_1_2[3], label="1/2, +1/2")
plot(ts, data_2_1[3], label="3/2, -3/2")
plot(ts, data_2_2[3], label="3/2, -1/2")
plot(ts, data_2_3[3], label="3/2, +1/2")
plot(ts, data_2_4[3], label="3/2, +3/2")
ylim([0, 1])
grid()
xlabel("t")
ylabel("p(S)")
legend(fontsize=8)
# NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

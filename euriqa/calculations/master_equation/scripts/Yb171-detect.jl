#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot
using LsqFit

include("../utils.jl")

const Γe = 2π * 19.6

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

    couple_F0 = true

    θ = π / 2

    Ωπ = Ω * cos(θ)
    Ωσ = Ω * sin(θ) / sqrt(2)
    # π transition
    # (0, 0) -> (1, 0)
    H_add_Ω!(H, Ωπ, 1, 7)
    # (1, -1) -> (1, -1)
    H_add_Ω!(H, -Ωπ, 2, 6)
    if couple_F0
        # (1, 0) -> (0, 0)
        H_add_Ω!(H, Ωπ, 3, 5)
    end
    # (1, 1) -> (1, 1)
    H_add_Ω!(H, Ωπ, 4, 8)

    # σ⁺ transition
    # (0, 0) -> (1, 1)
    H_add_Ω!(H, Ωσ, 1, 8)
    if couple_F0
        # (1, -1) -> (0, 0)
        H_add_Ω!(H, Ωσ, 2, 5)
    end
    # (1, -1) -> (1, 0)
    H_add_Ω!(H, Ωσ, 2, 7)
    # (1, 0) -> (1, 1)
    H_add_Ω!(H, -Ωσ, 3, 8)

    # σ⁻ transition
    # (0, 0) -> (1, -1)
    H_add_Ω!(H, Ωσ, 1, 6)
    # (1, 0) -> (1, -1)
    H_add_Ω!(H, Ωσ, 3, 6)
    if couple_F0
        # (1, 1) -> (0, 0)
        H_add_Ω!(H, Ωσ, 4, 5)
    end
    # (1, 1) -> (1, 0)
    H_add_Ω!(H, -Ωσ, 4, 7)
    add_H!(sys, H)

    # π transition
    Cπ = new_C(sys)
    # (0, 0) -> (1, 0)
    C_add_decay!(Cπ, Γe / 3, 1, 7)
    # (1, -1) -> (1, -1)
    C_add_decay!(Cπ, Γe / 3, 2, 6)
    # (1, 0) -> (0, 0)
    C_add_decay!(Cπ, Γe / 3, 3, 5)
    # (1, 1) -> (1, 1)
    C_add_decay!(Cπ, Γe / 3, 4, 8)
    add_C!(sys, Cπ)

    # σ⁺ transition
    Cσ⁺ = new_C(sys)
    # (0, 0) -> (1, 1)
    C_add_decay!(Cσ⁺, Γe / 3, 1, 8)
    # (1, -1) -> (0, 0)
    C_add_decay!(Cσ⁺, Γe / 3, 2, 5)
    # (1, -1) -> (1, 0)
    C_add_decay!(Cσ⁺, Γe / 3, 2, 7)
    # (1, 0) -> (1, 1)
    C_add_decay!(Cσ⁺, Γe / 3, 3, 8)
    add_C!(sys, Cσ⁺)

    # σ⁻ transition
    Cσ⁻ = new_C(sys)
    # (0, 0) -> (1, -1)
    C_add_decay!(Cσ⁻, Γe / 3, 1, 6)
    # (1, 0) -> (1, -1)
    C_add_decay!(Cσ⁻, Γe / 3, 3, 6)
    # (1, 1) -> (0, 0)
    C_add_decay!(Cσ⁻, Γe / 3, 4, 5)
    # (1, 1) -> (1, 0)
    C_add_decay!(Cσ⁻, Γe / 3, 4, 7)
    add_C!(sys, Cσ⁻)

    return sys
end

const ρ0 = zeros(8, 8)
ρ0[1, 1] = 1
const ρ1 = zeros(8, 8)
ρ1[3, 3] = 1

# model1(x, p) = p[1] .* (1 .- exp.(.-x .* (p[2] / p[1])))

# function compute_detection_data(s, Δ, B)
#     Ω = sqrt(s / 2) * Γe

#     sys = get_sys(Ω, Δ, B)
#     ts = Float64[]
#     p01s = Float64[]
#     p10s = Float64[]

#     t = 1.0
#     while true
#         ρ = propagate(sys, ρ1, t)
#         p10 = real(ρ[1, 1])
#         ρ = propagate(sys, ρ0, t)
#         p01 = 1 - real(ρ[1, 1])

#         push!(ts, t)
#         push!(p10s, p10)
#         push!(p01s, p01 .* 25)

#         if length(ts) > 5 && p10 > 0.4 && p01 > 0.02
#             break
#         end
#         t *= 2
#     end
#     fit10 = curve_fit(model1, ts, p10s, [1.0, 1 / t], lower=[0.0, 0.0])
#     r10 = fit10.param[2]
#     fit01 = curve_fit(model1, ts, p01s, [1.0, 1 / t], lower=[0.0, 0.0])
#     r01 = fit01.param[2] ./ 25

#     ts = range(0.5, 9.5, 20)
#     pe = 0.0
#     for t in ts
#         ρ = propagate(sys, ρ1, t)
#         pe += real(ρ[5, 5]) + real(ρ[6, 6]) + real(ρ[7, 7]) + real(ρ[8, 8])
#     end
#     γ = pe / length(ts) * Γe
#     return γ, r10, r01
# end

# function compute_detection_power_data(ss, Δ, B)
#     γs = Float64[]
#     r10s = Float64[]
#     r01s = Float64[]
#     for s in ss
#         γ, r10, r01 = compute_detection_data(s, Δ, B)
#         push!(γs, γ)
#         push!(r10s, r10)
#         push!(r01s, r01)
#     end
#     return γs, r10s, r01s
# end

const B0 = 3
# const B1 = 20.0

# const ss = range(0.01, 3, 100)

# const Bs = [2, 2.8, 3, 4, 5]
# const power_data = compute_detection_power_data.(Ref(ss), 0, Bs)

# figure()
# for (B, pd) in zip(Bs, power_data)
#     plot(ss, pd[1], label="B=$(B) G")
# end
# grid()
# xlim([0, ss[end]])
# xlabel("s")
# ylim([0, ylim()[2]])
# ylabel("Scattering Rate (MHz)")
# legend(fontsize=12, ncols=2)

# figure()
# for (B, pd) in zip(Bs, power_data)
#     plot(ss, pd[2] .* 1000, label="B=$(B) G")
# end
# grid()
# xlim([0, ss[end]])
# xlabel("s")
# ylim([0, ylim()[2]])
# ylabel("1\$\\rightarrow\$0 Rate (kHz)")
# legend(fontsize=12, ncols=2)

# figure()
# for (B, pd) in zip(Bs, power_data)
#     plot(ss, pd[3] .* 1000, label="B=$(B) G")
# end
# grid()
# xlim([0, ss[end]])
# xlabel("s")
# ylim([0, ylim()[2]])
# ylabel("0\$\\rightarrow\$1 Rate (kHz)")
# legend(fontsize=12, ncols=2)

# figure()
# for (B, pd) in zip(Bs, power_data)
#     plot(ss, sqrt.(pd[2] .* pd[3]) ./ pd[1] .* 100, label="B=$(B) G")
# end
# grid()
# xlim([0, ss[end]])
# xlabel("s")
# ylim([0, ylim()[2]])
# ylabel("\$\\sqrt{R_{0\\rightarrow1}R_{1\\rightarrow0}}/\\gamma\$ (%)")
# legend(fontsize=12, ncols=2)

const Δ0 = 2π * 2.105e3 * 2

const sys_1 = get_sys(2π * 1, Δ0, B0)
const sys_2 = get_sys(2π * 2, Δ0, B0)
const sys_5 = get_sys(2π * 5, Δ0, B0)
const sys_10 = get_sys(2π * 10, Δ0, B0)
const sys_20 = get_sys(2π * 20, Δ0, B0)
const sys_50 = get_sys(2π * 50, Δ0, B0)
const sys_100 = get_sys(2π * 100, Δ0, B0)
const sys_200 = get_sys(2π * 200, Δ0, B0)
const sys_500 = get_sys(2π * 500, Δ0, B0)

const ts = range(0, 30000, 1001)

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

const data_1 = @time get_plot_data(sys_1, ρ1, ts)
const data_2 = @time get_plot_data(sys_2, ρ1, ts)
const data_5 = @time get_plot_data(sys_5, ρ1, ts)
const data_10 = @time get_plot_data(sys_10, ρ1, ts)
const data_20 = @time get_plot_data(sys_20, ρ1, ts)
const data_50 = @time get_plot_data(sys_50, ρ1, ts)
const data_100 = @time get_plot_data(sys_100, ρ1, ts)
const data_200 = @time get_plot_data(sys_200, ρ1, ts)
const data_500 = @time get_plot_data(sys_500, ρ1, ts)

const prefix = joinpath(@__DIR__, "../imgs/three-states-scatter")

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

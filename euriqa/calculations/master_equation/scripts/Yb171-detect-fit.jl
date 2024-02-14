#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using NaCsData.Fitting
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

const ρ1 = zeros(8, 8)
ρ1[3, 3] = 1

function get_scatter_rate(s, Δ, B)
    Ω = sqrt(s / 2) * Γe

    sys = get_sys(Ω, Δ, B)

    ts = range(0.5, 149.5, 20)
    pe = 0.0
    for t in ts
        ρ = propagate(sys, ρ1, t)
        pe += real(ρ[5, 5]) + real(ρ[6, 6]) + real(ρ[7, 7]) + real(ρ[8, 8])
    end
    γ = pe / length(ts) * Γe
    return γ
end

function get_model(Δ, B)
    function model(x, p)
        return p[1] .* get_scatter_rate.(x.^2 .* p[2], Δ, B)
    end
end

const model = get_model(0, 2.8)
const amps = [0.32, 0.3, 0.28, 0.25, 0.2, 0.17, 0.15, 0.13, 0.1, 0.08]
const rates = [10.30, 11.14, 11.32, 11.38, 11.09, 9.32, 7.27, 6.41, 3.46, 2.07]

const fit = fit_data(model, amps, rates, [1.0, 10.0], plot_npts=200)

figure()
plot(amps, rates, "C0o")
plot(fit.plotx, fit.ploty, "C0")
grid()
xlim([0, 0.33])
ylim([0, 12.3])
xlabel("Amplitude")
ylabel("Count")
# NaCsPlot.maybe_save()

NaCsPlot.maybe_show()

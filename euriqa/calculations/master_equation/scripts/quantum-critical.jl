#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsPlot
using PyPlot

include("../utils.jl")

function get_sys(N, ω, Ω, λ, κ)
    sys = System{ComplexF64}(N * 2)

    H = new_H(sys)
    for i in 1:N
        n = i - 1
        H_add_E!(H, ω * n - Ω / 2, i)
        H_add_E!(H, ω * n + Ω / 2, i + N)
        if i >= N
            continue
        end
        H_add_Ω!(H, -λ * √(n + 1), i, i + 1 + N)
        H_add_Ω!(H, -λ * √(n + 1), i + 1, i + N)
    end
    add_H!(sys, H)

    C = new_C(sys)
    for i in 1:N - 1
        C_add_decay!(C, κ * √(i), i, i + 1)
        C_add_decay!(C, κ * √(i), i + N, i + 1 + N)
    end
    add_C!(sys, C)

    return sys
end

function propagate_n̄(sys, ρ0, t)
    ρ = @time propagate(sys, ρ0, t)
    N = size(ρ, 1) ÷ 2
    n̄ = 0.0
    for i in 1:N
        n = i - 1
        n̄ += n * real(ρ[i, i] + ρ[i + N, i + N])
    end
    return n̄
end

function propagate_n̄_λ0(λ, t)
    N = 16
    sys = get_sys(N, 0.1, 20.0, λ, 0.05)
    ρ0 = zeros(2N, 2N)
    ρ0[1, 1] = 1.0
    return propagate_n̄(sys, ρ0, t)
end

const λs = range(0, 3, 101)

# const N = 16
# const sys_0 = get_sys(N, 1.0, 1.0, 0.1, 1.0)
# const sys_1 = get_sys(N, 1.0, 1.0, 3.0, 1.0)

# const ts = range(0, 30, 51)
# const ρ0 = zeros(2N, 2N)
# ρ0[2, 2] = 1.0

const prefix = joinpath(@__DIR__, "../imgs/quantum-critical")

figure()
plot(λs, propagate_n̄_λ0.(λs, 300))
# plot(ts, (t->(sys_0, ρ0, t)).(ts), "C0")
# plot(ts, (t->propagate_n̄(sys_1, ρ0, t)).(ts), "C0")
# ylim([0, 1])
grid()
xlabel("λ")
ylabel("n")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

#!/usr/bin/julia

include("utils.jl")

const N_A = 6.02214076e23
const m_Yb171 = 170.9363258e-3 / N_A
const ω_m = 2π * 25e3
const λ = 578e-9

const η1 = Trap.η(m_Yb171, ω_m / (2π), 2π / λ)

@show η1

const Ω = 2π * 50e3

const T = 10 * 20e3 * 2π

function evolve(H, t, nbar)
    N = size(H.H, 1) ÷ 2
    pe = 0.0
    pg = 0.0
    nmax = ceil(Int, max(nbar * 4, 3))
    a = (1 + nbar) / nbar
    U = get_U(H, t)
    for n in 0:nmax
        w = a^(-n - 1) * (a - 1)
        pe′ = get_pe(U * get_ψ(N, n))
        pe += w * pe′
        pg += w * (1 - pe′)
    end
    return pe, pg
end

const N = 60
const H_drive = ExpCache(get_H(N, ω_m, 0, Ω, η1))

const ts = range(0, 20e-6, 1000)

const pes = Vector{Float64}(undef, length(ts))
const pgs = Vector{Float64}(undef, length(ts))

@time for i in 1:length(ts)
    pes[i], pgs[i] = evolve(H_drive, ts[i], T / ω_m)
end

figure()
plot(ts .* 1e6, pes, label="e")
plot(ts .* 1e6, pgs, label="g")
grid()
xlabel("t (\$\\mu s\$)")

NaCsPlot.maybe_show()

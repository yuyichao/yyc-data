#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: RandSetBits, RandDepol, Rand2QDepol
using NaCsPlot
using PyPlot

const Clf = NaCsCalc.Clifford

struct RNGs{T,RM,RD,RD2}
    rm::RM
    rd::RD
    rd2::RD2
    function RNGs{T}(p, pmeasure) where T
        rm = RandSetBits{T}(pmeasure)
        rd = RandDepol{T}(p)
        rd2 = Rand2QDepol{T}(p)
        return new{T,typeof(rm),typeof(rd),typeof(rd2)}(rm, rd, rd2)
    end
end

const stab_xs = ((1, 2, 3, 4), (2, 3, 5, 6), (3, 4, 6, 7))
const stab_zs = ((1, 2, 3, 4), (2, 3, 5, 6), (3, 4, 6, 7))

const logic_x = (5, 6, 7)
const logic_z = (5, 6, 7)

function init!(state::Clf.PauliString)
    empty!(state)
end

function apply_noisy_cnot!(state, i, j, rng)
    @inbounds Clf.apply!(state, Clf.CNOTGate(), i, j)
    x1, z1, x2, z2 = rand(rng.rd2)
    Clf.inject_pauli!(state, x1, z1, i)
    Clf.inject_pauli!(state, x2, z2, j)
end

function apply_transverse_cnot!(state, offset1, offset2, rng)
    for i in 1:7
        apply_noisy_cnot!(state, (offset1 - 1) * 7 + i,
                          (offset2 - 1) * 7 + i, rng)
    end
end

function inject_error!(state, rng)
    @inbounds for i in 1:7
        xmask, zmask = rand(rng.rd)
        Clf.inject_pauli!(state, xmask, zmask, i)
    end
end

@inline function measure_error(val, rng)
    return val ⊻ rand(rng.rm)
end

function decode_error!(state, stab_x_vals, stab_z_vals)
    T = eltype(state)
    @inbounds begin
        Clf.inject_pauli!(state, stab_z_vals[1] & ~stab_z_vals[2] & ~stab_z_vals[3],
                          zero(T), 1)
        Clf.inject_pauli!(state, ~stab_z_vals[1] & stab_z_vals[2] & ~stab_z_vals[3],
                          zero(T), 5)
        Clf.inject_pauli!(state, ~stab_z_vals[1] & ~stab_z_vals[2] & stab_z_vals[3],
                          zero(T), 7)

        Clf.inject_pauli!(state, stab_z_vals[1] & stab_z_vals[2] & ~stab_z_vals[3],
                          zero(T), 2)
        Clf.inject_pauli!(state, stab_z_vals[1] & ~stab_z_vals[2] & stab_z_vals[3],
                          zero(T), 4)
        Clf.inject_pauli!(state, ~stab_z_vals[1] & stab_z_vals[2] & stab_z_vals[3],
                          zero(T), 6)

        Clf.inject_pauli!(state, stab_z_vals[1] & stab_z_vals[2] & stab_z_vals[3],
                          zero(T), 3)

        Clf.inject_pauli!(state, zero(T),
                          stab_x_vals[1] & ~stab_x_vals[2] & ~stab_x_vals[3], 1)
        Clf.inject_pauli!(state, zero(T),
                          ~stab_x_vals[1] & stab_x_vals[2] & ~stab_x_vals[3], 5)
        Clf.inject_pauli!(state, zero(T),
                          ~stab_x_vals[1] & ~stab_x_vals[2] & stab_x_vals[3], 7)

        Clf.inject_pauli!(state, zero(T),
                          stab_x_vals[1] & stab_x_vals[2] & ~stab_x_vals[3], 2)
        Clf.inject_pauli!(state, zero(T),
                          stab_x_vals[1] & ~stab_x_vals[2] & stab_x_vals[3], 4)
        Clf.inject_pauli!(state, zero(T),
                          ~stab_x_vals[1] & stab_x_vals[2] & stab_x_vals[3], 6)

        Clf.inject_pauli!(state, zero(T),
                          stab_x_vals[1] & stab_x_vals[2] & stab_x_vals[3], 3)
    end
end

function correction_round!(state::Clf.PauliString, rng)
    T = eltype(state)
    @inbounds for i in 8:21
        state.xs[i] = zero(T)
        state.zs[i] = zero(T)
    end
    @inbounds for i in 8:21
        xmask, zmask = rand(rng.rd)
        Clf.inject_pauli!(state, xmask, zmask, i)
    end
    apply_transverse_cnot!(state, 1, 2, rng)
    apply_transverse_cnot!(state, 3, 1, rng)

    stab_x_vals = map(stab->@inbounds(measure_error(Clf.measure_stabilizer_x(state, stab .+ 14), rng)), stab_xs)
    stab_z_vals = map(stab->@inbounds(measure_error(Clf.measure_stabilizer_z(state, stab .+ 7), rng)), stab_zs)

    decode_error!(state, stab_x_vals, stab_z_vals)

    return
end

function correct_error!(state, rng)
    stab_x_vals = map(stab->@inbounds(Clf.measure_stabilizer_x(state, stab)), stab_xs)
    stab_z_vals = map(stab->@inbounds(Clf.measure_stabilizer_z(state, stab)), stab_zs)

    decode_error!(state, stab_x_vals, stab_z_vals)

    return
end

mutable struct ErrorStat
    err::Int
    tot::Int
    ErrorStat() = new(0, 0)
end

function collect_results(state::Clf.PauliString, err_stat::ErrorStat)
    err_stat.tot += Clf.nbits(eltype(state))
    err_stat.err += @inbounds Clf.count_ones(
        Clf.measure_stabilizer_z(state, logic_z) |
            Clf.measure_stabilizer_x(state, logic_x))
    return
end

function simulate(state, p, pmeasure, round, n)
    err_stat = ErrorStat()
    rng = RNGs{eltype(state)}(p, pmeasure)
    while err_stat.tot[] < n / sqrt(p)
        init!(state)
        inject_error!(state, rng)
        for i in 1:round
            correction_round!(state, rng)
        end
        correct_error!(state, rng)
        collect_results(state, err_stat)
    end
    return err_stat
end

function calc_errors(state, ps, pmeasure_scale, round, n)
    eps = Vector{Float64}(undef, length(ps))
    for (i, p) in enumerate(ps)
        err_stat = simulate(state, p, p * pmeasure_scale, round, n)
        eps[i] = err_stat.err / err_stat.tot
    end
    return eps
end

const ps = range(0.0001, 1, 200)
const ps2 = exp10.(range(-5, -1, 501))

const prefix = joinpath(@__DIR__, "imgs/steane_circuit/")

const state2 = Clf.PauliString{UInt128}(21)

const eps1 = (
    ps = ps,
    n0 = @time(calc_errors(state2, ps, 0, 0, 1280 * 256)),
    n1 = @time(calc_errors(state2, ps, 0, 1, 1280 * 256)),
    m1 = @time(calc_errors(state2, ps, 1, 1, 1280 * 256)),
    n4 = @time(calc_errors(state2, ps, 0, 4, 1280 * 256)),
    m4 = @time(calc_errors(state2, ps, 1, 4, 1280 * 256)),
    n16 = @time(calc_errors(state2, ps, 0, 16, 1280 * 256)),
    m16 = @time(calc_errors(state2, ps, 1, 16, 1280 * 256))
)

const eps2 = (
    ps = ps2,
    n0 = @time(calc_errors(state2, ps2, 0, 0, 320 * 256 * 8)),
    n1 = @time(calc_errors(state2, ps2, 0, 1, 320 * 256 * 8)),
    m1 = @time(calc_errors(state2, ps2, 1, 1, 320 * 256 * 8)),
    n4 = @time(calc_errors(state2, ps2, 0, 4, 320 * 256 * 3)),
    m4 = @time(calc_errors(state2, ps2, 1, 4, 320 * 256 * 3)),
    n16 = @time(calc_errors(state2, ps2, 0, 16, 320 * 256)),
    m16 = @time(calc_errors(state2, ps2, 1, 16, 320 * 256))
)

function plot_lines(eps)
    plot(eps.ps, eps.n0, "C1", label="0 round")
    plot(eps.ps, eps.n1, "C3", label="1 round")
    plot(eps.ps, eps.m1, "C4", label="1 round (measure error)")
    plot(eps.ps, eps.n4, "C5", label="4 rounds")
    plot(eps.ps, eps.m4, "C6", label="4 rounds (measure error)")
    plot(eps.ps, eps.n16, "C7", label="16 rounds")
    plot(eps.ps, eps.m16, "C8", label="16 rounds (measure error)")
    plot(eps.ps, eps.ps, "C0--")
end

function plot_lines_scaled(eps)
    plot(eps.ps, eps.n0, "C1", label="0 round")
    plot(eps.ps, eps.n1, "C3", label="1 round")
    plot(eps.ps, eps.m1, "C4", label="1 round (measure error)")
    plot(eps.ps .* 4, eps.n4, "C5", label="4 rounds")
    plot(eps.ps .* 4, eps.m4, "C6", label="4 rounds (measure error)")
    plot(eps.ps .* 16, eps.n16, "C7", label="16 rounds")
    plot(eps.ps .* 16, eps.m16, "C8", label="16 rounds (measure error)")
    plot(eps.ps, eps.ps, "C0--")
end

figure()
plot_lines(eps1)
legend(fontsize=12)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 1])
ylim([0, 1])
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)full")

figure()
plot_lines(eps2)
legend(fontsize=12)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 0.09])
ylim([0, 0.09])
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)zoomin")

figure()
plot_lines(eps2)
legend(fontsize=12)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([1e-5, 0.1])
ylim([1e-7, 0.11])
xscale("log")
yscale("log")
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)log")

figure()
plot_lines_scaled(eps2)
legend(fontsize=12)
grid()
xlabel("Physical Error (scaled)")
ylabel("Logical Error")
xlim([1e-5, 0.1])
ylim([1e-7, 0.11])
xscale("log")
yscale("log")
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)log_scaled")

NaCsPlot.maybe_show()

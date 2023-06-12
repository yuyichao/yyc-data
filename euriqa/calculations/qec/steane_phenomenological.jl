#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: rand_setbits
using NaCsPlot
using PyPlot

const Clf = NaCsCalc.Clifford

const stab_xs = ((1, 2, 3, 4), (2, 3, 5, 6), (3, 4, 6, 7))
const stab_zs = ((1, 2, 3, 4), (2, 3, 5, 6), (3, 4, 6, 7))

const logic_x = (5, 6, 7)
const logic_z = (5, 6, 7)

function init!(state::Clf.PauliString)
    empty!(state)
end

function inject_error!(state, p)
    Tele = eltype(state)
    @inbounds for i in 1:7
        xmask, zmask = Clf.rand_depol_error(Tele, p)
        Clf.inject_pauli!(state, xmask, zmask, i)
    end
end

function repeated_error(p, d)
    perror = 0.0
    count = 1.0
    for i in 0:((d - 1) ÷ 2)
        perror += (1 - p)^i * p^(d - i) * count
        count = count * (d - i) / (i + 1)
    end
    return min(perror, 1.0)
end

function measure_error(val::T, pmeasure) where T
    return val ⊻ rand_setbits(T, pmeasure)
end

function correct_error!(state, pmeasure)
    T = eltype(state)
    stab_x_vals = map(stab->@inbounds(measure_error(Clf.measure_stabilizer_x(state, stab), pmeasure)), stab_xs)
    stab_z_vals = map(stab->@inbounds(measure_error(Clf.measure_stabilizer_z(state, stab), pmeasure)), stab_zs)

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

function simulate_idle(state, p, d, n, pmeasure)
    err_stat = ErrorStat()
    while err_stat.tot[] < n
        init!(state)
        inject_error!(state, p)
        correct_error!(state, pmeasure)
        collect_results(state, err_stat)
    end
    return err_stat
end

function calc_errors(state, ps, d, n; pmeasure=nothing)
    eps = Vector{Float64}(undef, length(ps))
    for (i, p) in enumerate(ps)
        err_stat = simulate_idle(state, p, d, n,
                                 pmeasure === nothing ?
                                     repeated_error(p, d) : pmeasure)
        eps[i] = err_stat.err / err_stat.tot
    end
    return eps
end

const ps = [range(0, 0.1, 501); range(0.1, 1, 501)]

const prefix = joinpath(@__DIR__, "imgs/steane_phenomenological")

const state2 = Clf.PauliString{UInt128}(7)

const eps0 = @time calc_errors(state2, ps, 1, 6400 * 256, pmeasure=0)
const eps1 = @time calc_errors(state2, ps, 1, 6400 * 256)
const eps3 = @time calc_errors(state2, ps, 3, 6400 * 256)
const eps5 = @time calc_errors(state2, ps, 5, 6400 * 256)
const eps7 = @time calc_errors(state2, ps, 7, 6400 * 256)
const eps101 = @time calc_errors(state2, ps, 101, 6400 * 256)

function plot_lines()
    plot(ps, eps1, "C2", label="measure 1 round")
    plot(ps, eps3, "C3", label="measure 3 rounds")
    plot(ps, eps5, "C4", label="measure 5 rounds")
    plot(ps, eps7, "C5", label="measure 7 rounds")
    plot(ps, eps101, "C6", label="measure 101 rounds")
    plot(ps, eps0, "C1", alpha=0.4, label="no measure error")
    plot(ps, ps, "C0--")
end

figure()
plot_lines()
legend(fontsize=12)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 1])
ylim([0, 1])
title("[[7, 1, 3]] Phenomenological")
NaCsPlot.maybe_save(prefix)

figure()
plot_lines()
legend(fontsize=12)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0.02, 0.12])
ylim([0.007, 0.13])
title("[[7, 1, 3]] Phenomenological")
NaCsPlot.maybe_save("$(prefix)_zoomin")

figure()
plot_lines()
legend(fontsize=12)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0.02, 0.3])
ylim([0.007, 0.8])
xscale("log")
yscale("log")
title("[[7, 1, 3]] Phenomenological")
NaCsPlot.maybe_save("$(prefix)_log")

NaCsPlot.maybe_show()

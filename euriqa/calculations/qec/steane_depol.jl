#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: RandDepol
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

function inject_error!(state, rd)
    Tele = eltype(state)
    @inbounds for i in 1:7
        xmask, zmask = rand(rd)
        Clf.inject_pauli!(state, xmask, zmask, i)
    end
end

function correct_error!(state)
    T = eltype(state)
    stab_x_vals = map(stab->@inbounds(Clf.measure_stabilizer_x(state, stab)), stab_xs)
    stab_z_vals = map(stab->@inbounds(Clf.measure_stabilizer_z(state, stab)), stab_zs)

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

function simulate_idle(state, p, n)
    err_stat = ErrorStat()
    rd = RandDepol{eltype(state)}(p)
    while err_stat.tot[] < n
        init!(state)
        inject_error!(state, rd)
        correct_error!(state)
        collect_results(state, err_stat)
    end
    return err_stat
end

function calc_errors(state, ps, n)
    eps = Vector{Float64}(undef, length(ps))
    for (i, p) in enumerate(ps)
        err_stat = simulate_idle(state, p, n)
        eps[i] = err_stat.err / err_stat.tot
    end
    return eps
end

const ps = range(0, 1, 501)

const prefix = joinpath(@__DIR__, "imgs/steane_depol")

const state1 = Clf.StabilizerState(7)
const state2 = Clf.PauliString{UInt128}(7)

const eps_diff = @time calc_errors(state2, ps, 6400 * 256)

figure()
plot(ps, eps_diff, "C2")
plot(ps, ps, "C0--")
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 1])
ylim([0, 1])
title("Steane code [[7, 1, 3]] Depolarizing")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

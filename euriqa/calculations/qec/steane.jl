#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: RandSetBits
using NaCsPlot
using PyPlot

const Clf = NaCsCalc.Clifford

const stab_xs = ((1, 2, 3, 4), (2, 3, 5, 6), (3, 4, 6, 7))
const stab_zs = ((1, 2, 3, 4), (2, 3, 5, 6), (3, 4, 6, 7))

const logic_x = (5, 6, 7)
const logic_z = (5, 6, 7)

function init!(state::Clf.StabilizerState, isx, v)
    if isx
        Clf.init_state_x!(state)
    else
        Clf.init_state_z!(state)
    end
    for stab in stab_xs
        v2, det = Clf.measure_xs!(state, stab; force=false)
        @assert det == isx
        @assert !v2
    end
    for stab in stab_zs
        v2, det = Clf.measure_zs!(state, stab; force=false)
        @assert det != isx
        @assert !v2
    end
    if v
        if isx
            for z in logic_z
                Clf.apply!(state, Clf.ZGate(), z)
            end
        else
            for x in logic_x
                Clf.apply!(state, Clf.XGate(), x)
            end
        end
    end
end

function init!(state::Clf.PauliString, isx, v)
    empty!(state)
end

function inject_error!(state, sb)
    @inbounds for i in 1:7
        Clf.inject_pauli!(state, rand(sb), rand(sb), i)
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
    xerr::Int
    xtot::Int
    zerr::Int
    ztot::Int
    ErrorStat() = new(0, 0, 0, 0)
end

function collect_result_x(state, v, err_stat::ErrorStat)
    err_stat.xtot += Clf.nbits(eltype(state))
    err_stat.xerr += @inbounds Clf.count_ones(
        Clf.measure_stabilizer_x(state, logic_x, v))
    return
end

function collect_result_z(state, v, err_stat::ErrorStat)
    err_stat.ztot += Clf.nbits(eltype(state))
    err_stat.zerr += @inbounds Clf.count_ones(
        Clf.measure_stabilizer_z(state, logic_z, v))
    return
end

function collect_results(state::Clf.StabilizerState, isx, v, err_stat::ErrorStat)
    if isx
        collect_result_x(state, v, err_stat)
    else
        collect_result_z(state, v, err_stat)
    end
end

function collect_results(state::Clf.PauliString, isx, v, err_stat::ErrorStat)
    collect_result_x(state, false, err_stat)
    collect_result_z(state, false, err_stat)
end

function simulate_idle(state, p, n)
    err_stat = ErrorStat()
    sb = RandSetBits{eltype(state)}(p)
    while err_stat.xtot[] + err_stat.ztot[] < 4 * n
        for isx in (false, true)
            for v in (false, true)
                init!(state, isx, v)
                inject_error!(state, sb)
                correct_error!(state)
                collect_results(state, isx, v, err_stat)
            end
        end
    end
    return err_stat
end

function calc_errors(state, ps, n)
    xps = Vector{Float64}(undef, length(ps))
    zps = Vector{Float64}(undef, length(ps))
    for (i, p) in enumerate(ps)
        err_stat = simulate_idle(state, p, n)
        xps[i] = err_stat.xerr / err_stat.xtot
        zps[i] = err_stat.zerr / err_stat.ztot
    end
    return xps, zps
end

const ps = range(0, 1, 501)

const prefix = joinpath(@__DIR__, "imgs/steane")

const state1 = Clf.StabilizerState(7)
const state2 = Clf.PauliString{UInt128}(7)

const xps_state, zps_state = @time calc_errors(state1, ps, 3200)
const xps_diff, zps_diff = @time calc_errors(state2, ps, 3200 * 256)

figure()
plot(ps, xps_state, "C0", label="x (state)", alpha=0.3)
plot(ps, xps_diff, "C2", label="x (error)")
plot(ps, zps_state, "C1", label="z (state)", alpha=0.3)
plot(ps, zps_diff, "C3", label="z (error)")
plot(ps, ps, "C0--")
legend(fontsize=13, ncol=2)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 1])
ylim([0, 1])
title("Steane code [[7, 1, 3]]")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

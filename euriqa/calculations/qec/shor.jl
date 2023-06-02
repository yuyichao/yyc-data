#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: rand_setbits
using NaCsPlot
using PyPlot

const Clf = NaCsCalc.Clifford

const stab_xs = ((1, 2), (2, 3), (4, 5), (5, 6), (7, 8), (8, 9))
const stab_zs = ((1, 2, 3, 4, 5, 6), (4, 5, 6, 7, 8, 9))

const logic_x = (1, 4, 7)
const logic_z = (1, 2, 3)

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

function inject_error!(state, px, pz)
    Tele = eltype(state)
    for i in 1:9
        Clf.inject_pauli!(state, rand_setbits(Tele, px), rand_setbits(Tele, pz), i)
    end
end

function correct_error!(state)
    T = eltype(state)
    stab_x_vals = map(stab->Clf.measure_stabilizer_x(state, stab), stab_xs)
    stab_z_vals = map(stab->Clf.measure_stabilizer_z(state, stab), stab_zs)

    Clf.inject_pauli!(state, stab_z_vals[1] & ~stab_z_vals[2], zero(T), 1)
    Clf.inject_pauli!(state, stab_z_vals[1] & stab_z_vals[2], zero(T), 4)
    Clf.inject_pauli!(state, ~stab_z_vals[1] & stab_z_vals[2], zero(T), 7)

    Clf.inject_pauli!(state, zero(T), stab_x_vals[1] & ~stab_x_vals[2], 1)
    Clf.inject_pauli!(state, zero(T), stab_x_vals[1] & stab_x_vals[2], 2)
    Clf.inject_pauli!(state, zero(T), ~stab_x_vals[1] & stab_x_vals[2], 3)

    Clf.inject_pauli!(state, zero(T), stab_x_vals[3] & ~stab_x_vals[4], 4)
    Clf.inject_pauli!(state, zero(T), stab_x_vals[3] & stab_x_vals[4], 5)
    Clf.inject_pauli!(state, zero(T), ~stab_x_vals[3] & stab_x_vals[4], 6)

    Clf.inject_pauli!(state, zero(T), stab_x_vals[5] & ~stab_x_vals[6], 7)
    Clf.inject_pauli!(state, zero(T), stab_x_vals[5] & stab_x_vals[6], 8)
    Clf.inject_pauli!(state, zero(T), ~stab_x_vals[5] & stab_x_vals[6], 9)

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
    err_stat.xerr += Clf.count_ones(Clf.measure_stabilizer_x(state, logic_x, v))
    return
end

function collect_result_z(state, v, err_stat::ErrorStat)
    err_stat.ztot += Clf.nbits(eltype(state))
    err_stat.zerr += Clf.count_ones(Clf.measure_stabilizer_z(state, logic_z, v))
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

function simulate_idle(state, px, pz, n)
    err_stat = ErrorStat()
    for isx in (false, true)
        for v in (false, true)
            for i in 1:n
                init!(state, isx, v)
                inject_error!(state, px, pz)
                correct_error!(state)
                collect_results(state, isx, v, err_stat)
            end
        end
    end
    return err_stat
end

function calc_errors(state, ps, n)
    xps = Float64[]
    zps = Float64[]
    for p in ps
        err_stat = simulate_idle(state, p, p, n)
        push!(xps, err_stat.xerr / err_stat.xtot)
        push!(zps, err_stat.zerr / err_stat.ztot)
    end
    return xps, zps
end

const ps = range(0, 1, 501)

const prefix = joinpath(@__DIR__, "imgs/shor")

const state1 = Clf.StabilizerState(9)
const state2 = Clf.PauliString{UInt}(9)

const xps_state, zps_state = @time calc_errors(state1, ps, 3200 * 2)
const xps_diff, zps_diff = @time calc_errors(state2, ps, 3200)

figure()
plot(ps, xps_state, label="x (state)")
plot(ps, zps_state, label="z (state)")
plot(ps, xps_diff, label="x (error)")
plot(ps, zps_diff, label="z (error)")
plot(ps, ps, "C0--")
legend(fontsize=13, ncol=2)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 1])
ylim([0, 1])
title("Shor code [[9, 1, 3]]")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

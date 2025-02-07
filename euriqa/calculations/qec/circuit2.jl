#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: RandSetBits, RandDepol, Rand2QDepol
using NaCsPlot
using PyPlot

const Clf = NaCsCalc.Clifford

struct ErrorBit
    bit::Int
    type::Int
end

function _compute_syndrome(errs::Vector{ErrorBit}, stab_x, stab_z)
    res = false
    for err in errs
        if err.type == 1
            res ⊻= stab_z[err.bit]
        elseif err.type == 2
            res ⊻= stab_x[err.bit] ⊻ stab_z[err.bit]
        elseif err.type == 3
            res ⊻= stab_x[err.bit]
        end
    end
    return res
end

function _next_err_type!(errs)
    for i in 1:length(errs)
        err = errs[i]
        if err.type >= 3
            errs[i] = ErrorBit(err.bit, 1)
            continue
        end
        errs[i] = ErrorBit(err.bit, err.type + 1)
        return true
    end
    return false
end

function _next_bit!(errs, nq)
    prev_bit = 0
    nerr = length(errs)
    for i in 1:nerr
        err = errs[nerr - i + 1]
        max_bit = nq + 1 - i
        if err.bit < max_bit
            bit_base = err.bit
            for j in 1:i
                errs[nerr - i + j] = ErrorBit(bit_base + j, 1)
            end
            return true
        end
    end
    return false
end

function decoding_table(nq, stabs_x, stabs_z)
    nstab = length(stabs_x)
    lot = Dict{Vector{Bool},Vector{ErrorBit}}()
    @assert length(stabs_z) == nstab
    for nerr in 1:nq
        errs = [ErrorBit(i, 1) for i in 1:nerr]
        while true
            syndrome = [_compute_syndrome(errs, stabs_x[i], stabs_z[i])
                        for i in 1:nstab]
            if any(syndrome) && !(syndrome in keys(lot))
                lot[syndrome] = copy(errs)
                if length(lot) == 2^nstab - 1
                    break
                end
            end
            if !_next_err_type!(errs)
                if !_next_bit!(errs, nq)
                    break
                end
            end
        end
    end
    return lot
end

function init!(state::Clf.PauliString)
    empty!(state)
    return
end

function noise_1q!(state, i, rd)
    xmask, zmask = rand(rd)
    @inbounds Clf.inject_pauli!(state, xmask, zmask, i)
    return
end

function noise_2q!(state, i, j, rd2)
    x1, z1, x2, z2 = rand(rd2)
    @inbounds Clf.inject_pauli!(state, x1, z1, i)
    @inbounds Clf.inject_pauli!(state, x2, z2, j)
    return
end

function apply_noisy!(state, gate::Clf.Clifford1Q, i, rd)
    @inbounds Clf.apply!(state, gate, i)
    noise_1q!(state, i, rd)
    return
end

function apply_noisy!(state, gate::Clf.Clifford2Q, i, j, rd2)
    @inbounds Clf.apply!(state, gate, i, j)
    noise_2q!(state, i, j, rd2)
    return
end

measure_noisy_z(state, i, rm) = Clf.measure_stabilizer_z(state, (i,)) ⊻ rand(rm)

function correct_error!(state, lot, stab_vals)
    T = eltype(state)
    for (syndrome, errs) in lot
        mask = ~zero(T)
        for (s, sv) in zip(syndrome, stab_vals)
            if s
                mask &= sv
            else
                mask &= ~sv
            end
        end
        if mask != 0
            for err in errs
                if err.type == 1
                    Clf.inject_pauli!(state, mask, zero(T), err.bit)
                elseif err.type == 2
                    Clf.inject_pauli!(state, mask, mask, err.bit)
                elseif err.type == 3
                    Clf.inject_pauli!(state, zero(T), mask, err.bit)
                end
            end
        end
    end
    return
end

struct UniformInit{T,RD}
    rd::RD
    function UniformInit{T}(p) where T
        rd = RandDepol{T}(p)
        return new{T,typeof(rd)}(rd)
    end
end

function init!(state, init::UniformInit)
    init!(state)
    for i in 1:state.n
        noise_1q!(state, i, init.rd)
    end
end

struct UniformInit2{T,RD}
    n1::Int
    rd::RD
    rd_2::RD
    function UniformInit2{T}(n1, p1, p2) where T
        rd = RandDepol{T}(p1)
        rd_2 = RandDepol{T}(p2)
        return new{T,typeof(rd)}(n1, rd, rd_2)
    end
end

function init!(state, init::UniformInit2)
    init!(state)
    for i in 1:init.n1
        noise_1q!(state, i, init.rd)
    end
    for i in (init.n1 + 1):state.n
        noise_1q!(state, i, init.rd_2)
    end
end

mutable struct ErrorStat
    err::Int
    tot::Int
    ErrorStat() = new(0, 0)
end

function collect_results(state::Clf.PauliString, err_stat::ErrorStat,
                         logics_x, logics_z)
    err_stat.tot += Clf.nbits(eltype(state))
    mask = zero(eltype(state))
    for (logic_x, logic_z) in zip(logics_x, logics_z)
        mask |= @inbounds Clf.measure_stabilizer(state, logic_x, logic_z)
    end
    err_stat.err += @inbounds Clf.count_ones(mask)
    return
end

struct RawStabMeasureCircuit{T,Init,RM,RD,RD2}
    nq::Int
    init::Init
    stabs_x::Vector{Vector{Bool}}
    stabs_z::Vector{Vector{Bool}}
    logics_x::Vector{Vector{Bool}}
    logics_z::Vector{Vector{Bool}}
    rms::Vector{RM}
    rds::Vector{RD}
    rd2s::Matrix{RD2}
    function RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
                                      pms, pds, pd2s, init::Init) where {T,Init}
        nstab = length(stabs_x)
        @assert nstab > 0
        @assert length(stabs_z) == nstab
        nq = length(stabs_x[1])
        @assert all(length.(stabs_x) .== nq)
        @assert all(length.(stabs_z) .== nq)
        nl = length(logics_x)
        @assert length(logics_z) == nl
        @assert all(length.(logics_x) .== nq)
        @assert all(length.(logics_z) .== nq)

        @assert length(pms) == nstab
        @assert length(pds) == nstab
        @assert size(pd2s) == (nq, nstab)

        rms = [RandSetBits{T}(p) for p in pms]
        rds = [RandDepol{T}(p) for p in pds]
        rd2s = [Rand2QDepol{T}(p) for p in pd2s]
        return new{T,Init,eltype(rms),eltype(rds),eltype(rd2s)}(
            nq, init, stabs_x, stabs_z, logics_x, logics_z, rms, rds, rd2s)
    end
end

function _measure_stab_raw(circ::RawStabMeasureCircuit, state, stabi)
    anc = stabi + circ.nq
    rd = circ.rds[stabi]
    stab_x = circ.stabs_x[stabi]
    stab_z = circ.stabs_z[stabi]
    apply_noisy!(state, Clf.HGate(), anc, rd)
    for i in 1:circ.nq
        x = stab_x[i]
        z = stab_z[i]
        rd2 = circ.rd2s[i, stabi]
        if x
            if z
                apply_noisy!(state, Clf.CYGate(), anc, i, rd2)
            else
                apply_noisy!(state, Clf.CXGate(), anc, i, rd2)
            end
        elseif z
            apply_noisy!(state, Clf.CZGate(), anc, i, rd2)
        end
    end
    apply_noisy!(state, Clf.HGate(), anc, rd)
    return measure_noisy_z(state, anc, circ.rms[stabi])
end

function run(circ::RawStabMeasureCircuit{T}, nrep, final_round=true) where T
    nstab = length(circ.stabs_x)
    state = Clf.PauliString{T}(circ.nq + nstab)
    err_stat = ErrorStat()
    stab_vals = zeros(T, nstab)
    lot = decoding_table(circ.nq, circ.stabs_x, circ.stabs_z)
    while err_stat.tot[] < nrep
        init!(state, circ.init)
        stab_vals .= _measure_stab_raw.(Ref(circ), Ref(state), 1:nstab)
        correct_error!(state, lot, stab_vals)
        if final_round
            stab_vals .= Clf.measure_stabilizer.(Ref(state), stabs_x, stabs_z)
            correct_error!(state, lot, stab_vals)
        end
        collect_results(state, err_stat, circ.logics_x, circ.logics_z)
    end
    return err_stat
end

function is_css_code(stabs_x, stabs_z)
    for (stab_x, stab_z) in zip(stabs_x, stabs_z)
        if any(stab_x) && any(stab_z)
            return false
        end
    end
    return true
end

struct SteaneMeasureCircuit{T,Init,RM,RD,RD2}
    nq::Int
    init::Init
    stabs_x::Vector{Vector{Bool}}
    stabs_z::Vector{Vector{Bool}}
    logics_x::Vector{Vector{Bool}}
    logics_z::Vector{Vector{Bool}}
    rms::Matrix{RM}
    rds::Matrix{RD}
    rd2s::Matrix{RD2}
    function SteaneMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
                                      pms, pds, pd2s, init::Init) where {T,Init}
        @assert is_css_code(stabs_x, stabs_z)

        nstab = length(stabs_x)
        @assert nstab > 0
        @assert length(stabs_z) == nstab
        nq = length(stabs_x[1])
        @assert all(length.(stabs_x) .== nq)
        @assert all(length.(stabs_z) .== nq)
        nl = length(logics_x)
        @assert length(logics_z) == nl
        @assert all(length.(logics_x) .== nq)
        @assert all(length.(logics_z) .== nq)

        @assert size(pms) == (nq, 2)
        @assert size(pds) == (nq, 2)
        @assert size(pd2s) == (nq, 2)

        rms = [RandSetBits{T}(p) for p in pms]
        rds = [RandDepol{T}(p) for p in pds]
        rd2s = [Rand2QDepol{T}(p) for p in pd2s]
        return new{T,Init,eltype(rms),eltype(rds),eltype(rd2s)}(
            nq, init, stabs_x, stabs_z, logics_x, logics_z, rms, rds, rd2s)
    end
end

function _measure_steane(circ::SteaneMeasureCircuit, state, copy_id, isz, res)
    offset = copy_id * circ.nq
    for i in 1:circ.nq
        if isz
            apply_noisy!(state, Clf.CNOTGate(), i, offset + i, circ.rd2s[i, copy_id])
        else
            rd = circ.rds[i, copy_id]
            apply_noisy!(state, Clf.HGate(), offset + i, rd)
            apply_noisy!(state, Clf.CNOTGate(), offset + i, i, circ.rd2s[i, copy_id])
            apply_noisy!(state, Clf.HGate(), offset + i, rd)
        end
        res[i] = measure_noisy_z(state, offset + i, circ.rms[i, copy_id])
    end
end

function _measure_stab_steane(circ::SteaneMeasureCircuit{T},
                              state, stabi, res_x, res_z) where T
    stab_x = circ.stabs_x[stabi]
    stab_z = circ.stabs_z[stabi]
    res = zero(T)
    for i in 1:circ.nq
        x = stab_x[i]
        z = stab_z[i]
        if x
            res ⊻= res_x[i]
        elseif z
            res ⊻= res_z[i]
        end
    end
    return res
end

function run(circ::SteaneMeasureCircuit{T}, nrep, final_round=true) where T
    nstab = length(circ.stabs_x)
    state = Clf.PauliString{T}(circ.nq * 3)
    err_stat = ErrorStat()
    res_x = zeros(T, circ.nq)
    res_z = zeros(T, circ.nq)
    stab_vals = zeros(T, nstab)
    lot = decoding_table(circ.nq, circ.stabs_x, circ.stabs_z)
    while err_stat.tot[] < nrep
        init!(state, circ.init)
        _measure_steane(circ, state, 1, false, res_x)
        _measure_steane(circ, state, 2, true, res_z)
        stab_vals .= _measure_stab_steane.(Ref(circ), Ref(state), 1:nstab,
                                           Ref(res_x), Ref(res_z))
        correct_error!(state, lot, stab_vals)
        if final_round
            stab_vals .= Clf.measure_stabilizer.(Ref(state), stabs_x, stabs_z)
            correct_error!(state, lot, stab_vals)
        end
        collect_results(state, err_stat, circ.logics_x, circ.logics_z)
    end
    return err_stat
end

# Steane
# const stabs_x = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false]]
# const stabs_z = [[false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]

# const stabs_x = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]
# const stabs_z = [[false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]

# const stabs_x = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false]]
# const stabs_z = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]

# const logics_x = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, false, false, false]]
# const logics_z = [[false, false, false, false, false, false, false],
#                   [false, false, false, false, true, true, true]]

# const logics_x = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, true, true, true]]
# const logics_z = [[false, false, false, false, false, false, false],
#                   [false, false, false, false, true, true, true]]

# const logics_x = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, false, false, false]]
# const logics_z = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, true, true, true]]

const stabs_x = [[true, true, false, false, false, false, false, false, false],
                 [false, true, true, false, false, false, false, false, false],
                 [false, false, false, true, true, false, false, false, false],
                 [false, false, false, false, true, true, false, false, false],
                 [false, false, false, false, false, false, true, true, false],
                 [false, false, false, false, false, false, false, true, true],
                 [false, false, false, false, false, false, false, false, false],
                 [false, false, false, false, false, false, false, false, false]]
const stabs_z = [[false, false, false, false, false, false, false, false, false],
                 [false, false, false, false, false, false, false, false, false],
                 [false, false, false, false, false, false, false, false, false],
                 [false, false, false, false, false, false, false, false, false],
                 [false, false, false, false, false, false, false, false, false],
                 [false, false, false, false, false, false, false, false, false],
                 [true, true, true, true, true, true, false, false, false],
                 [false, false, false, true, true, true, true, true, true]]

const logics_x = [[true, false, false, true, false, false, true, false, false],
                  [false, false, false, false, false, false, false, false, false]]
const logics_z = [[false, false, false, false, false, false, false, false, false],
                  [true, true, true, false, false, false, false, false, false]]


function calc_errors(ps, n)
    eps = Vector{Float64}(undef, length(ps))
    T = UInt128
    nstab = length(stabs_x)
    nq = length(stabs_x[1])
    for (i, p) in enumerate(ps)
        # init = UniformInit2{T}(nq, p, 0.0)
        init = UniformInit{T}(p)
        circ = RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
                                        zeros(nstab), zeros(nstab),
                                        zeros(nq, nstab), init)
        # circ = RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
        #                                 ones(nstab) .* 0.0001, ones(nstab) .* 0.0001,
        #                                 ones(nq, nstab) .* 0.0001, init)
        # circ = RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
        #                                 ones(nstab) .* p, ones(nstab) .* p,
        #                                 ones(nq, nstab) .* p, init)
        # circ = SteaneMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
        #                                zeros(nq, 2), zeros(nq, 2),
        #                                zeros(nq, 2), init)
        # circ = SteaneMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
        #                                ones(nq, 2) .* 0.0001, ones(nq, 2) .* 0.0001,
        #                                ones(nq, 2) .* 0.0001, init)
        # circ = SteaneMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
        #                                ones(nq, 2) .* p, ones(nq, 2) .* p,
        #                                ones(nq, 2) .* p, init)
        err_stat = run(circ, n)
        eps[i] = err_stat.err / err_stat.tot
    end
    return eps
end

const ps = range(0.0001, 1, 200)
const ps2 = exp10.(range(-5, 0, 1001))

const prefix = joinpath(@__DIR__, "imgs/steane_circuit/")

const eps2 = (
    ps = ps2,
    n0 = @time(calc_errors(ps2, 320 * 256 * 8)),
)

function plot_lines(eps)
    plot(eps.ps, eps.n0, "C1")
    plot(eps.ps, eps.ps, "C0--")
end

figure()
plot_lines(eps2)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 0.15])
ylim([0, 0.15])
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)zoomin")

figure()
plot_lines(eps2)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([1e-5, 0.15])
ylim([1e-7, 0.11])
xscale("log")
yscale("log")
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)log")

NaCsPlot.maybe_show()

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

function apply_noisy!(state, gate::Clf.Clifford1Q, i, rd::RandDepol)
    @inbounds Clf.apply!(state, gate, i)
    noise_1q!(state, i, rd)
    return
end

function apply_noisy!(state::Clf.PauliString{T}, gate::Clf.Clifford1Q,
                      i, err::NTuple{2,T}) where T
    @inbounds Clf.apply!(state, gate, i)
    @inbounds Clf.inject_pauli!(state, err..., i)
    return
end

function apply_noisy!(state, gate::Clf.Clifford2Q, i, j, rd2::Rand2QDepol)
    @inbounds Clf.apply!(state, gate, i, j)
    noise_2q!(state, i, j, rd2)
    return
end

function apply_noisy!(state::Clf.PauliString{T}, gate::Clf.Clifford2Q,
                      i, j, err::NTuple{4,T}) where T
    @inbounds Clf.apply!(state, gate, i, j)
    @inbounds Clf.inject_pauli!(state, err[1], err[2], i)
    @inbounds Clf.inject_pauli!(state, err[3], err[4], j)
    return
end

measure_noisy_z(state, i, rm::RandSetBits) =
    Clf.measure_stabilizer_z(state, (i,)) ⊻ rand(rm)

measure_noisy_z(state::Clf.PauliString{T}, i, err::T) where T =
    Clf.measure_stabilizer_z(state, (i,)) ⊻ err

mutable struct RNGs{T,RM,RD,RD2}
    const rms::Vector{RM}
    const rds::Vector{RD}
    const rd2s::Vector{RD2}
    idx_1q::Int
    idx_2q::Int
    idx_m::Int
    function RNGs{T}(pms, pds, pd2s) where {T}
        rms = [RandSetBits{T}(p) for p in pms]
        rds = [RandDepol{T}(p) for p in pds]
        rd2s = [Rand2QDepol{T}(p) for p in pd2s]
        return new{T,eltype(rms),eltype(rds),eltype(rd2s)}(rms, rds, rd2s, 1, 1, 1)
    end
end

function init!(rngs::RNGs)
    rngs.idx_1q = 1
    rngs.idx_2q = 1
    rngs.idx_m = 1
    return
end

function apply_noisy_next!(state, gate::Clf.Clifford1Q, i, rngs)
    idx = rngs.idx_1q
    rngs.idx_1q = idx + 1
    apply_noisy!(state, gate, i, rngs.rds[idx])
    return
end

function apply_noisy_next!(state, gate::Clf.Clifford2Q, i, j, rngs)
    idx = rngs.idx_2q
    rngs.idx_2q = idx + 1
    apply_noisy!(state, gate, i, j, rngs.rd2s[idx])
    return
end

function measure_noisy_z_next(state, i, rngs)
    idx = rngs.idx_m
    rngs.idx_m = idx + 1
    return measure_noisy_z(state, i, rngs.rms[idx])
end

mutable struct Errors{T}
    const ems::Vector{T}
    const eds::Vector{NTuple{2,T}}
    const ed2s::Vector{NTuple{4,T}}
    idx_1q::Int
    idx_2q::Int
    idx_m::Int
    function Errors{T}(ems, eds, ed2s) where {T}
        return new{T}(ems, eds, ed2s, 1, 1, 1)
    end
end

function init!(errs::Errors)
    errs.idx_1q = 1
    errs.idx_2q = 1
    errs.idx_m = 1
    return
end

function apply_noisy_next!(state, gate::Clf.Clifford1Q, i, errs)
    idx = errs.idx_1q
    errs.idx_1q = idx + 1
    apply_noisy!(state, gate, i, errs.rds[idx])
    return
end

function apply_noisy_next!(state, gate::Clf.Clifford2Q, i, j, errs)
    idx = errs.idx_2q
    errs.idx_2q = idx + 1
    apply_noisy!(state, gate, i, j, errs.rd2s[idx])
    return
end

function measure_noisy_z_next(state, i, errs)
    idx = errs.idx_m
    errs.idx_m = idx + 1
    return measure_noisy_z(state, i, errs.rms[idx])
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

function _convert_stab_orders(nq, stabs_x, stabs_z, stab_orders::Nothing)
    nstab = length(stabs_x)
    return [(i, j) for i in 1:nq, j in 1:nstab if stabs_x[j][i] || stabs_z[j][i]]
end

function _check_stab_orders(nq, stabs_x, stabs_z, stab_orders)
    nstab = length(stabs_x)
    masks = zeros(Bool, nq, nstab)
    for (q, s) in stab_orders
        masks[q, s] = true
    end
    for q in 1:nq
        for s in 1:nstab
            @assert masks[q, s] == stabs_x[s][q] || stabs_z[s][q]
        end
    end
end

function _convert_stab_orders(nq, stabs_x, stabs_z,
                              stab_orders::AbstractVector{VI}
                              where VI <: AbstractVector{Int})
    nstab = length(stabs_x)
    new_orders = NTuple{2,Int}[]
    for s in 1:nstab
        for q in stab_orders[s]
            push!(new_orders, (q, s))
        end
    end
    _check_stab_orders(nq, stabs_x, stabs_z, new_orders)
    return new_orders
end

function _convert_stab_orders(nq, stabs_x, stabs_z,
                              stab_orders::AbstractVector{NTuple{2,Int}})
    _check_stab_orders(nq, stabs_x, stabs_z, stab_orders)
    return stab_orders
end

struct RawStabMeasureCircuit{T,Init,R}
    nq::Int
    init::Init
    stabs_x::Vector{Vector{Bool}}
    stabs_z::Vector{Vector{Bool}}
    logics_x::Vector{Vector{Bool}}
    logics_z::Vector{Vector{Bool}}
    stab_orders::Vector{NTuple{2,Int}}
    rngs::R
    function RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
                                      rngs::R, init::Init;
                                      stab_orders=nothing) where {T,Init,R}
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

        stab_orders = _convert_stab_orders(nq, stabs_x, stabs_z, stab_orders)

        return new{T,Init,R}(nq, init, stabs_x, stabs_z, logics_x, logics_z,
                             stab_orders, rngs)
    end
end

function _do_stab_raw_2q!(circ::RawStabMeasureCircuit, state)
    for (i, stabi) in circ.stab_orders
        anc = stabi + circ.nq
        stab_x = circ.stabs_x[stabi]
        stab_z = circ.stabs_z[stabi]
        x = stab_x[i]
        z = stab_z[i]
        if x
            if z
                apply_noisy_next!(state, Clf.CYGate(), anc, i, circ.rngs)
            else
                apply_noisy_next!(state, Clf.CXGate(), anc, i, circ.rngs)
            end
        elseif z
            apply_noisy_next!(state, Clf.CZGate(), anc, i, circ.rngs)
        end
    end
end

function run(circ::RawStabMeasureCircuit{T}, nrep, final_round=true) where T
    nstab = length(circ.stabs_x)
    state = Clf.PauliString{T}(circ.nq + nstab)
    err_stat = ErrorStat()
    stab_vals = zeros(T, nstab)
    lot = decoding_table(circ.nq, circ.stabs_x, circ.stabs_z)
    while err_stat.tot[] < nrep
        init!(state, circ.init)
        init!(circ.rngs)
        for stabi in 1:nstab
            anc = stabi + circ.nq
            apply_noisy_next!(state, Clf.HGate(), anc, circ.rngs)
        end
        _do_stab_raw_2q!(circ, state)
        for stabi in 1:nstab
            anc = stabi + circ.nq
            apply_noisy_next!(state, Clf.HGate(), anc, circ.rngs)
            stab_vals[stabi] = measure_noisy_z_next(state, anc, circ.rngs)
        end

        correct_error!(state, lot, stab_vals)
        if final_round
            stab_vals .= Clf.measure_stabilizer.(Ref(state),
                                                 circ.stabs_x, circ.stabs_z)
            correct_error!(state, lot, stab_vals)
        end
        collect_results(state, err_stat, circ.logics_x, circ.logics_z)
    end
    return err_stat
end

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

function decoding_table(nq, t, stabs_x, stabs_z)
    nstab = length(stabs_x)
    lot = Dict{Vector{Bool},Vector{ErrorBit}}()
    @assert length(stabs_z) == nstab
    for nerr in 1:t
        errs = [ErrorBit(i, 1) for i in 1:nerr]
        while true
            syndrome = [_compute_syndrome(errs, stabs_x[i], stabs_z[i])
                        for i in 1:nstab]
            if any(syndrome) && !(syndrome in keys(lot))
                lot[syndrome] = copy(errs)
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

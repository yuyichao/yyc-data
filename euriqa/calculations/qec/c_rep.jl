#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc
using NaCsCalc.Utils: rand_setbits
using NaCsPlot
using PyPlot

struct State{T}
    vals::Vector{T}
    State{T}(n) where T = new(zeros(T, n))
end

function init!(state::State{T}, v) where T
    state.vals .= v
    return state
end

function inject_error!(state::State{T}, p) where T
    @inbounds for i in 1:length(state.vals)
        state.vals[i] ⊻= rand_setbits(T, p)
    end
end

function decode(state::State{T}) where T
    nbits = T === Bool ? 1 : sizeof(T) * 8
    counters = ntuple(i->Ref(0), Val(nbits))
    threshold = length(state.vals) ÷ 2
    for i in 0:(nbits - 1)
        c = 0
        for v in state.vals
            c += (v >> i) & 1
        end
        counters[i + 1][] = c
    end
    if T === Bool
        return counters[1][] > threshold
    end
    res = zero(T)
    for i in 0:(nbits - 1)
        if counters[i + 1][] > threshold
            res |= one(T) << i
        end
    end
    return res
end

function simulate_idle(state::State{T}, p, n) where T
    nbits = T === Bool ? 1 : sizeof(T) * 8
    ntotal = 0
    nerror = 0
    for i in 1:n
        v = rand(T)
        init!(state, v)
        inject_error!(state, p)
        v2 = decode(state)
        err = v ⊻ v2
        ntotal += nbits
        nerror += T === Bool ? err : count_ones(err)
    end
    return nerror, ntotal
end

function calc_error(n, p)
    state = State{UInt}(n)
    e, t = simulate_idle(state, p, 10000 ÷ sizeof(UInt))
    return e / t
end

const ps = range(0, 1, 1001)

const prefix = joinpath(@__DIR__, "imgs/c_rep")

figure()
plot(ps, calc_error.(3, ps), label="3")
plot(ps, calc_error.(5, ps), label="5")
plot(ps, calc_error.(9, ps), label="9")
plot(ps, calc_error.(11, ps), label="11")
plot(ps, calc_error.(13, ps), label="13")
plot(ps, calc_error.(15, ps), label="15")
plot(ps, calc_error.(17, ps), label="17")
plot(ps, ps, "C0--")
legend(fontsize=13, ncol=2)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 1])
ylim([0, 1])
title("Classical repetition code")
NaCsPlot.maybe_save(prefix)

NaCsPlot.maybe_show()

#!/usr/bin/julia

using NaCsCalc
using BenchmarkTools

const Clf = NaCsCalc.Clifford

function test_bacon_shor(n)
    data = Matrix{Int}(undef, n, n)
    nqubit = 0

    for i in 1:n
        for j in 1:n
            nqubit += 1
            data[i, j] = nqubit
        end
    end
    state = Clf.StabilizerState(nqubit)
    # logical |+> preparation
    @inbounds for i in 1:n
        d1 = data[i, 1]
        Clf.apply!(state, Clf.HGate(), d1)
        for j in 1:n
            Clf.apply!(state, Clf.CNOTGate(), d1, data[i, j])
        end
    end

    # println(state)

    # X gauge measurement
    @inbounds for i in 1:n - 1
        parity = false
        for j in 1:n
            meas, = Clf.measure_xs!(state, (data[i, j], data[i + 1, j]))
            parity = parity ⊻ meas
        end
        @assert !parity
    end

    # Z gauge measurement
    @inbounds for j in 1:n - 1
        parity = false
        for i in 1:n
            meas, = Clf.measure_zs!(state, (data[i, j], data[i, j + 1]))
            parity = parity ⊻ meas
        end
        @assert !parity
    end
end

# 2.297 ms (7 allocations: 91.12 KiB)
@btime test_bacon_shor(20)
# 166.944 ms (8 allocations: 1.26 MiB)
@btime test_bacon_shor(40)
# 1.251 s (9 allocations: 6.29 MiB)
@btime test_bacon_shor(60)
# 5.886 s (9 allocations: 19.68 MiB)
@btime test_bacon_shor(80)
# 39.800 s (9 allocations: 47.99 MiB)
@btime test_bacon_shor(100)

#!/usr/bin/julia

using NaCsCalc
using BenchmarkTools

const Clf = NaCsCalc.Clifford

function test_bacon_shor(::Type{SST}, n) where SST
    data = Matrix{Int}(undef, n, n)
    nqubit = 0

    for i in 1:n
        for j in 1:n
            nqubit += 1
            data[i, j] = nqubit
        end
    end
    state = SST(nqubit)
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

# 1.492 ms (7 allocations: 91.12 KiB)
# 1.574 ms (5 allocations: 104.38 KiB)
@btime test_bacon_shor(Clf.StabilizerState, 20)
@btime test_bacon_shor(Clf.InvStabilizerState, 20)
# 76.296 ms (8 allocations: 1.26 MiB)
# 30.600 ms (5 allocations: 1.29 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 40)
@btime test_bacon_shor(Clf.InvStabilizerState, 40)
# 678.306 ms (9 allocations: 6.29 MiB)
# 162.199 ms (6 allocations: 6.41 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 60)
@btime test_bacon_shor(Clf.InvStabilizerState, 60)
# 2.919 s (9 allocations: 19.68 MiB)
# 651.168 ms (6 allocations: 19.59 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 80)
@btime test_bacon_shor(Clf.InvStabilizerState, 80)
# 9.585 s (9 allocations: 47.99 MiB)
# 2.276 s (7 allocations: 48.32 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 100)
@btime test_bacon_shor(Clf.InvStabilizerState, 100)

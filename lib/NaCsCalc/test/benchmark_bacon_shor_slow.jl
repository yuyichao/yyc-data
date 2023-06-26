#!/usr/bin/julia

using NaCsCalc
using BenchmarkTools

const Clf = NaCsCalc.Clifford

function test_bacon_shor(n)
    data = Matrix{Int}(undef, n, n)
    anc_z = Matrix{Int}(undef, n, n - 1)
    anc_x = Matrix{Int}(undef, n - 1, n)
    nqubit = 0

    for j in 1:n
        for i in 1:n
            nqubit += 1
            data[i, j] = nqubit
        end
    end
    for j in 1:(n - 1)
        for i in 1:n
            nqubit += 1
            anc_z[i, j] = nqubit
        end
    end
    for j in 1:n
        for i in 1:(n - 1)
            nqubit += 1
            anc_x[i, j] = nqubit
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
            anc = anc_x[i, j]
            Clf.apply!(state, Clf.HGate(), anc)
            Clf.apply!(state, Clf.CNOTGate(), anc, data[i, j])
            Clf.apply!(state, Clf.CNOTGate(), anc, data[i + 1, j])
            Clf.apply!(state, Clf.HGate(), anc)
            meas, = Clf.measure_z!(state, anc)
            parity = parity ⊻ meas
        end
        @assert !parity
    end

    # Z gauge measurement
    @inbounds for j in 1:n - 1
        parity = false
        for i in 1:n
            anc = anc_z[i, j]
            Clf.apply!(state, Clf.CNOTGate(), anc, data[i, j])
            Clf.apply!(state, Clf.CNOTGate(), anc, data[i, j + 1])
            meas, = Clf.measure_z!(state, anc)
            parity = parity ⊻ meas
        end
        @assert !parity
    end
end

# 19.133 ms (9 allocations: 683.02 KiB)
@btime test_bacon_shor(20)
# 847.656 ms (9 allocations: 10.71 MiB)
@btime test_bacon_shor(40)
# 22.424 s (13 allocations: 54.53 MiB)
@btime test_bacon_shor(60)
# 152.242 s (13 allocations: 173.05 MiB)
@btime test_bacon_shor(80)

#!/usr/bin/julia

using NaCsCalc
using BenchmarkTools

const Clf = NaCsCalc.Clifford

function test_bacon_shor(::Type{SST}, n) where SST
    data = Matrix{Int}(undef, n, n)
    anc_z = Matrix{Int}(undef, n, n - 1)
    anc_x = Matrix{Int}(undef, n - 1, n)
    nqubit = 0

    for i in 1:n
        for j in 1:n
            nqubit += 1
            data[i, j] = nqubit
        end
    end
    for i in 1:n
        for j in 1:(n - 1)
            nqubit += 1
            anc_z[i, j] = nqubit
        end
    end
    for i in 1:(n - 1)
        for j in 1:n
            nqubit += 1
            anc_x[i, j] = nqubit
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

# 5.549 ms (10 allocations: 698.75 KiB)
# 2.906 ms (7 allocations: 737.39 KiB)
@btime test_bacon_shor(Clf.StabilizerState, 20)
@btime test_bacon_shor(Clf.InvStabilizerState, 20)
# 96.414 ms (10 allocations: 10.77 MiB)
# 52.857 ms (7 allocations: 10.71 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 40)
@btime test_bacon_shor(Clf.InvStabilizerState, 40)
# 702.525 ms (13 allocations: 54.68 MiB)
# 344.018 ms (11 allocations: 54.86 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 60)
@btime test_bacon_shor(Clf.InvStabilizerState, 60)
# 3.669 s (13 allocations: 173.30 MiB)
# 2.309 s (11 allocations: 173.34 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 80)
@btime test_bacon_shor(Clf.InvStabilizerState, 80)
# 14.109 s (13 allocations: 424.48 MiB)
# 11.864 s (11 allocations: 424.08 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 100)
@btime test_bacon_shor(Clf.InvStabilizerState, 100)
# 55.807 s (13 allocations: 881.35 MiB)
# 36.468 s (11 allocations: 881.44 MiB)
@btime test_bacon_shor(Clf.StabilizerState, 120)
@btime test_bacon_shor(Clf.InvStabilizerState, 120)

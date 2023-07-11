#!/usr/bin/julia

using NaCsCalc
using BenchmarkTools

const Clf = NaCsCalc.Clifford

@inline function rand_2bits(n)
    v1 = rand(1:n)
    v2 = rand(1:(n - 1))
    if v2 >= v1
        v2 += 1
    end
    return v1, v2
end

function test_size_rep(state, rep)
    n = state.n
    Clf.init_state_z!(state)
    @inbounds for i in 1:rep
        Clf.measure_xs!(state, rand_2bits(n))
        Clf.measure_zs!(state, rand_2bits(n))
    end
end

function test_size(::Type{SST}, n, rep) where SST
    state = SST(n)
    @btime test_size_rep($state, $rep)
end

# 8.792 μs (0 allocations: 0 bytes)
# 12.083 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 4, 100)
test_size(Clf.InvStabilizerState, 4, 100)
# 22.416 μs (0 allocations: 0 bytes)
# 46.667 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 8, 200)
test_size(Clf.InvStabilizerState, 8, 200)
# 65.375 μs (0 allocations: 0 bytes)
# 173.666 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 16, 400)
test_size(Clf.InvStabilizerState, 16, 400)
# 212.083 μs (0 allocations: 0 bytes)
# 617.001 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 32, 800)
test_size(Clf.InvStabilizerState, 32, 800)
# 977.459 μs (0 allocations: 0 bytes)
# 2.339 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 64, 1600)
test_size(Clf.InvStabilizerState, 64, 1600)
# 12.327 ms (0 allocations: 0 bytes)
# 27.333 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 512, 1600)
test_size(Clf.InvStabilizerState, 512, 1600)
# 27.423 ms (0 allocations: 0 bytes)
# 71.413 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 512, 3200)
test_size(Clf.InvStabilizerState, 512, 3200)
# 568.799 ms (0 allocations: 0 bytes)
# 1.342 s (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 2048, 6400)
test_size(Clf.InvStabilizerState, 2048, 6400)

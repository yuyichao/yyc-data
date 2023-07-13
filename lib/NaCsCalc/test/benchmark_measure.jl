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

# 8.749 μs (0 allocations: 0 bytes)
# 9.875 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 4, 100)
test_size(Clf.InvStabilizerState, 4, 100)
# 22.541 μs (0 allocations: 0 bytes)
# 39.165 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 8, 200)
test_size(Clf.InvStabilizerState, 8, 200)
# 65.497 μs (0 allocations: 0 bytes)
# 139.619 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 16, 400)
test_size(Clf.InvStabilizerState, 16, 400)
# 213.074 μs (0 allocations: 0 bytes)
# 522.603 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 32, 800)
test_size(Clf.InvStabilizerState, 32, 800)
# 956.417 μs (0 allocations: 0 bytes)
# 1.985 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 64, 1600)
test_size(Clf.InvStabilizerState, 64, 1600)
# 12.282 ms (0 allocations: 0 bytes)
# 16.959 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 512, 1600)
test_size(Clf.InvStabilizerState, 512, 1600)
# 27.209 ms (0 allocations: 0 bytes)
# 45.283 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 512, 3200)
test_size(Clf.InvStabilizerState, 512, 3200)
# 616.076 ms (0 allocations: 0 bytes)
# 468.468 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 2048, 6400)
test_size(Clf.InvStabilizerState, 2048, 6400)

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

function test_size_rep(state, n, rep)
    for i in 1:rep
        Clf.measure_xs!(state, rand_2bits(n))
        Clf.measure_zs!(state, rand_2bits(n))
    end
end

function test_size(::Type{SST}, n, rep) where SST
    state = SST(n)
    @btime test_size_rep($state, $n, $rep)
end

# 9.458 μs (4 allocations: 384 bytes)
test_size(Clf.StabilizerState, 4, 100)
test_size(Clf.InvStabilizerState, 4, 100)
# 24.665 μs (4 allocations: 512 bytes)
test_size(Clf.StabilizerState, 8, 200)
test_size(Clf.InvStabilizerState, 8, 200)
# 70.330 μs (4 allocations: 784 bytes)
test_size(Clf.StabilizerState, 16, 400)
test_size(Clf.InvStabilizerState, 16, 400)
# 227.823 μs (4 allocations: 1.28 KiB)
test_size(Clf.StabilizerState, 32, 800)
test_size(Clf.InvStabilizerState, 32, 800)
# 1.000 ms (4 allocations: 3.27 KiB)
test_size(Clf.StabilizerState, 64, 1600)
test_size(Clf.InvStabilizerState, 64, 1600)
# 12.429 ms (6 allocations: 136.41 KiB)
test_size(Clf.StabilizerState, 512, 1600)
test_size(Clf.InvStabilizerState, 512, 1600)
# 27.677 ms (6 allocations: 136.41 KiB)
test_size(Clf.StabilizerState, 512, 3200)
test_size(Clf.InvStabilizerState, 512, 3200)
# 580.491 ms (7 allocations: 2.03 MiB)
test_size(Clf.StabilizerState, 2048, 6400)
test_size(Clf.InvStabilizerState, 2048, 6400)

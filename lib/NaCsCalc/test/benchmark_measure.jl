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

function test_size(n, rep)
    state = Clf.StabilizerState(n)
    for i in 1:rep
        Clf.measure_xs!(state, rand_2bits(n))
        Clf.measure_zs!(state, rand_2bits(n))
    end
end

# 9.458 μs (4 allocations: 384 bytes)
@btime test_size(4, 100)
# 24.874 μs (4 allocations: 512 bytes)
@btime test_size(8, 200)
# 72.955 μs (4 allocations: 784 bytes)
@btime test_size(16, 400)
# 232.325 μs (4 allocations: 1.28 KiB)
@btime test_size(32, 800)
# 1.077 ms (4 allocations: 3.27 KiB)
@btime test_size(64, 1600)
# 19.939 ms (6 allocations: 136.41 KiB)
@btime test_size(512, 1600)
# 45.787 ms (6 allocations: 136.41 KiB)
@btime test_size(512, 3200)
# 1.140 s (7 allocations: 2.03 MiB)
@btime test_size(2048, 6400)

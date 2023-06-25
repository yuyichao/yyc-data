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

# 9.458 μs (4 allocations: 320 bytes)
@btime test_size(4, 100)
# 24.835 μs (4 allocations: 400 bytes)
@btime test_size(8, 200)
# 71.462 μs (4 allocations: 544 bytes)
@btime test_size(16, 400)
# 224.177 μs (4 allocations: 864 bytes)
@btime test_size(32, 800)
# 1.026 ms (4 allocations: 2.39 KiB)
@btime test_size(64, 1600)
# 19.196 ms (6 allocations: 129.34 KiB)
@btime test_size(512, 1600)
# 44.490 ms (6 allocations: 129.34 KiB)
@btime test_size(512, 3200)
# 1.152 s (6 allocations: 2.00 MiB)
@btime test_size(2048, 6400)

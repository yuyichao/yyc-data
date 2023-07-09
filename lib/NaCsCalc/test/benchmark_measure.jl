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
    for i in 1:rep
        Clf.measure_xs!(state, rand_2bits(n))
        Clf.measure_zs!(state, rand_2bits(n))
    end
end

function test_size(::Type{SST}, n, rep) where SST
    state = SST(n)
    @btime test_size_rep($state, $rep)
end

# 8.999 μs (0 allocations: 0 bytes)
# 17.083 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 4, 100)
test_size(Clf.InvStabilizerState, 4, 100)
# 23.499 μs (0 allocations: 0 bytes)
# 94.166 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 8, 200)
test_size(Clf.InvStabilizerState, 8, 200)
# 66.577 μs (0 allocations: 0 bytes)
# 407.796 μs (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 16, 400)
test_size(Clf.InvStabilizerState, 16, 400)
# 215.855 μs (0 allocations: 0 bytes)
# 1.539 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 32, 800)
test_size(Clf.InvStabilizerState, 32, 800)
# 961.165 μs (0 allocations: 0 bytes)
# 5.862 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 64, 1600)
test_size(Clf.InvStabilizerState, 64, 1600)
# 15.016 ms (0 allocations: 0 bytes)
# 61.126 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 512, 1600)
test_size(Clf.InvStabilizerState, 512, 1600)
# 30.375 ms (0 allocations: 0 bytes)
# 126.432 ms (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 512, 3200)
test_size(Clf.InvStabilizerState, 512, 3200)
# 947.757 ms (0 allocations: 0 bytes)
# 2.664 s (0 allocations: 0 bytes)
test_size(Clf.StabilizerState, 2048, 6400)
test_size(Clf.InvStabilizerState, 2048, 6400)

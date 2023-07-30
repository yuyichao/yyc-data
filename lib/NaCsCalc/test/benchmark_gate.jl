#!/usr/bin/julia

using NaCsCalc
using BenchmarkTools

const Clf = NaCsCalc.Clifford

@inline function benchmark(state, gate::Clf.Clifford1Q)
    n = state.n
    @inbounds for i in 1:n
        Clf.apply!(state, gate, i)
    end
end

@inline function benchmark(state, gate::Clf.Clifford2Q)
    n = state.n
    @inbounds for i in 1:n
        for j in 1:n - 1
            j = j < i ? j : j + 1
            Clf.apply!(state, gate, i, j)
        end
    end
end

const state1_1 = Clf.StabilizerState(29)
const state1_2 = Clf.StabilizerState(271)
const state1_3 = Clf.StabilizerState(5000)
const state2_1 = Clf.InvStabilizerState(29)
const state2_2 = Clf.InvStabilizerState(271)
const state2_3 = Clf.InvStabilizerState(5000)

function benchmark_all(gate)
    println(gate)
    @btime benchmark($(state1_1), $gate)
    @btime benchmark($(state1_2), $gate)
    @btime benchmark($(state1_3), $gate)
    @btime benchmark($(state2_1), $gate)
    @btime benchmark($(state2_2), $gate)
    @btime benchmark($(state2_3), $gate)
end

# 50.234 ns (0 allocations: 0 bytes)
# 879.120 ns (0 allocations: 0 bytes)
# 169.243 μs (0 allocations: 0 bytes)
# 42.758 ns (0 allocations: 0 bytes)
# 827.050 ns (0 allocations: 0 bytes)
# 172.910 μs (0 allocations: 0 bytes)
println("Init X")
@btime Clf.init_state_x!($(state1_1))
@btime Clf.init_state_x!($(state1_2))
@btime Clf.init_state_x!($(state1_3))
@btime Clf.init_state_x!($(state2_1))
@btime Clf.init_state_x!($(state2_2))
@btime Clf.init_state_x!($(state2_3))

# 50.234 ns (0 allocations: 0 bytes)
# 878.300 ns (0 allocations: 0 bytes)
# 170.618 μs (0 allocations: 0 bytes)
# 42.800 ns (0 allocations: 0 bytes)
# 850.203 ns (0 allocations: 0 bytes)
# 166.619 μs (0 allocations: 0 bytes)
println("Init Z")
@btime Clf.init_state_z!($(state1_1))
@btime Clf.init_state_z!($(state1_2))
@btime Clf.init_state_z!($(state1_3))
@btime Clf.init_state_z!($(state2_1))
@btime Clf.init_state_z!($(state2_2))
@btime Clf.init_state_z!($(state2_3))

# 2.291 ns (0 allocations: 0 bytes)
# 2.250 ns (0 allocations: 0 bytes)
# 2.250 ns (0 allocations: 0 bytes)
# 2.291 ns (0 allocations: 0 bytes)
# 2.625 ns (0 allocations: 0 bytes)
# 2.250 ns (0 allocations: 0 bytes)
benchmark_all(Clf.IGate())

# 63.476 ns (0 allocations: 0 bytes)
# 2.426 μs (0 allocations: 0 bytes)
# 441.567 μs (0 allocations: 0 bytes)
# 38.349 ns (0 allocations: 0 bytes)
# 845.800 ns (0 allocations: 0 bytes)
# 266.073 μs (0 allocations: 0 bytes)
benchmark_all(Clf.HGate())

# 63.603 ns (0 allocations: 0 bytes)
# 1.483 μs (0 allocations: 0 bytes)
# 129.995 μs (0 allocations: 0 bytes)
# 7.291 ns (0 allocations: 0 bytes)
# 10.301 ns (0 allocations: 0 bytes)
# 101.137 ns (0 allocations: 0 bytes)
benchmark_all(Clf.XGate())

# 63.772 ns (0 allocations: 0 bytes)
# 1.429 μs (0 allocations: 0 bytes)
# 250.740 μs (0 allocations: 0 bytes)
# 8.382 ns (0 allocations: 0 bytes)
# 14.904 ns (0 allocations: 0 bytes)
# 182.187 ns (0 allocations: 0 bytes)
benchmark_all(Clf.YGate())

# 63.560 ns (0 allocations: 0 bytes)
# 1.475 μs (0 allocations: 0 bytes)
# 129.786 μs (0 allocations: 0 bytes)
# 7.632 ns (0 allocations: 0 bytes)
# 10.301 ns (0 allocations: 0 bytes)
# 70.979 ns (0 allocations: 0 bytes)
benchmark_all(Clf.ZGate())

# 63.305 ns (0 allocations: 0 bytes)
# 2.333 μs (0 allocations: 0 bytes)
# 292.239 μs (0 allocations: 0 bytes)
# 52.483 ns (0 allocations: 0 bytes)
# 999.900 ns (0 allocations: 0 bytes)
# 256.033 μs (0 allocations: 0 bytes)
benchmark_all(Clf.SGate())

# 63.262 ns (0 allocations: 0 bytes)
# 2.338 μs (0 allocations: 0 bytes)
# 289.573 μs (0 allocations: 0 bytes)
# 53.805 ns (0 allocations: 0 bytes)
# 1.008 μs (0 allocations: 0 bytes)
# 256.365 μs (0 allocations: 0 bytes)
benchmark_all(Clf.ISGate())

# 63.264 ns (0 allocations: 0 bytes)
# 2.324 μs (0 allocations: 0 bytes)
# 298.412 μs (0 allocations: 0 bytes)
# 52.484 ns (0 allocations: 0 bytes)
# 979.143 ns (0 allocations: 0 bytes)
# 270.537 μs (0 allocations: 0 bytes)
benchmark_all(Clf.SXGate())

# 63.263 ns (0 allocations: 0 bytes)
# 2.338 μs (0 allocations: 0 bytes)
# 299.370 μs (0 allocations: 0 bytes)
# 53.806 ns (0 allocations: 0 bytes)
# 1.008 μs (0 allocations: 0 bytes)
# 260.079 μs (0 allocations: 0 bytes)
benchmark_all(Clf.ISXGate())

# 63.477 ns (0 allocations: 0 bytes)
# 2.440 μs (0 allocations: 0 bytes)
# 442.993 μs (0 allocations: 0 bytes)
# 38.683 ns (0 allocations: 0 bytes)
# 845.229 ns (0 allocations: 0 bytes)
# 269.412 μs (0 allocations: 0 bytes)
benchmark_all(Clf.SYGate())

# 63.434 ns (0 allocations: 0 bytes)
# 2.431 μs (0 allocations: 0 bytes)
# 437.617 μs (0 allocations: 0 bytes)
# 38.557 ns (0 allocations: 0 bytes)
# 845.229 ns (0 allocations: 0 bytes)
# 269.287 μs (0 allocations: 0 bytes)
benchmark_all(Clf.ISYGate())

# 1.871 μs (0 allocations: 0 bytes)
# 645.435 μs (0 allocations: 0 bytes)
# 2.248 s (0 allocations: 0 bytes)
# 3.015 μs (0 allocations: 0 bytes)
# 494.899 μs (0 allocations: 0 bytes)
# 2.001 s (0 allocations: 0 bytes)
benchmark_all(Clf.CNOTGate())

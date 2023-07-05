#!/usr/bin/julia

using NaCsCalc
using BenchmarkTools

const Clf = NaCsCalc.Clifford

@inline function benchmark(state, gate::Clf.Clifford1Q)
    n = state.n
    for i in 1:n
        Clf.apply!(state, gate, i)
    end
end

@inline function benchmark(state, gate::Clf.Clifford2Q)
    n = state.n
    for i in 1:n
        for j in 1:n - 1
            j = j < i ? j : j + 1
            Clf.apply!(state, gate, i, j)
        end
    end
end

const state1 = Clf.StabilizerState(29)
const state2 = Clf.StabilizerState(271)
const state3 = Clf.StabilizerState(5000)

# 2.250 ns (0 allocations: 0 bytes)
# 2.291 ns (0 allocations: 0 bytes)
# 2.250 ns (0 allocations: 0 bytes)
println("IGate")
@btime benchmark($(state1), Clf.IGate())
@btime benchmark($(state2), Clf.IGate())
@btime benchmark($(state3), Clf.IGate())

# 64.029 ns (0 allocations: 0 bytes)
# 2.361 μs (0 allocations: 0 bytes)
# 437.651 μs (0 allocations: 0 bytes)
println("HGate")
@btime benchmark($(state1), Clf.HGate())
@btime benchmark($(state2), Clf.HGate())
@btime benchmark($(state3), Clf.HGate())

# 63.815 ns (0 allocations: 0 bytes)
# 1.421 μs (0 allocations: 0 bytes)
# 135.662 μs (0 allocations: 0 bytes)
println("XGate")
@btime benchmark($(state1), Clf.XGate())
@btime benchmark($(state2), Clf.XGate())
@btime benchmark($(state3), Clf.XGate())

# 64.520 ns (0 allocations: 0 bytes)
# 1.375 μs (0 allocations: 0 bytes)
# 222.534 μs (0 allocations: 0 bytes)
println("YGate")
@btime benchmark($(state1), Clf.YGate())
@btime benchmark($(state2), Clf.YGate())
@btime benchmark($(state3), Clf.YGate())

# 63.773 ns (0 allocations: 0 bytes)
# 1.425 μs (0 allocations: 0 bytes)
# 136.620 μs (0 allocations: 0 bytes)
println("ZGate")
@btime benchmark($(state1), Clf.ZGate())
@btime benchmark($(state2), Clf.ZGate())
@btime benchmark($(state3), Clf.ZGate())

# 65.474 ns (0 allocations: 0 bytes)
# 2.116 μs (0 allocations: 0 bytes)
# 296.656 μs (0 allocations: 0 bytes)
println("SGate")
@btime benchmark($(state1), Clf.SGate())
@btime benchmark($(state2), Clf.SGate())
@btime benchmark($(state3), Clf.SGate())

# 65.510 ns (0 allocations: 0 bytes)
# 2.111 μs (0 allocations: 0 bytes)
# 291.698 μs (0 allocations: 0 bytes)
println("ISGate")
@btime benchmark($(state1), Clf.ISGate())
@btime benchmark($(state2), Clf.ISGate())
@btime benchmark($(state3), Clf.ISGate())

# 65.540 ns (0 allocations: 0 bytes)
# 2.120 μs (0 allocations: 0 bytes)
# 292.198 μs (0 allocations: 0 bytes)
println("SXGate")
@btime benchmark($(state1), Clf.SGate())
@btime benchmark($(state2), Clf.SGate())
@btime benchmark($(state3), Clf.SGate())

# 65.516 ns (0 allocations: 0 bytes)
# 2.093 μs (0 allocations: 0 bytes)
# 294.490 μs (0 allocations: 0 bytes)
println("ISXGate")
@btime benchmark($(state1), Clf.ISGate())
@btime benchmark($(state2), Clf.ISGate())
@btime benchmark($(state3), Clf.ISGate())

# 1.946 μs (0 allocations: 0 bytes)
# 548.731 μs (0 allocations: 0 bytes)
# 2.252 s (0 allocations: 0 bytes)
println("CNOTGate")
@btime benchmark($(state1), Clf.CNOTGate())
@btime benchmark($(state2), Clf.CNOTGate())
@btime benchmark($(state3), Clf.CNOTGate())

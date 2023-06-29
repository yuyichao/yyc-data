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

# 56.485 ns (0 allocations: 0 bytes)
# 460.858 ns (0 allocations: 0 bytes)
println("IGate")
@btime benchmark($(state1), Clf.IGate())
@btime benchmark($(state2), Clf.IGate())

# 102.659 ns (0 allocations: 0 bytes)
# 1.954 μs (0 allocations: 0 bytes)
println("HGate")
@btime benchmark($(state1), Clf.HGate())
@btime benchmark($(state2), Clf.HGate())

# 82.037 ns (0 allocations: 0 bytes)
# 1.433 μs (0 allocations: 0 bytes)
println("XGate")
@btime benchmark($(state1), Clf.XGate())
@btime benchmark($(state2), Clf.XGate())

# 94.228 ns (0 allocations: 0 bytes)
# 1.521 μs (0 allocations: 0 bytes)
println("YGate")
@btime benchmark($(state1), Clf.YGate())
@btime benchmark($(state2), Clf.YGate())

# 83.806 ns (0 allocations: 0 bytes)
# 1.396 μs (0 allocations: 0 bytes)
println("ZGate")
@btime benchmark($(state1), Clf.ZGate())
@btime benchmark($(state2), Clf.ZGate())

# 97.893 ns (0 allocations: 0 bytes)
# 1.968 μs (0 allocations: 0 bytes)
println("SGate")
@btime benchmark($(state1), Clf.SGate())
@btime benchmark($(state2), Clf.SGate())

# 97.205 ns (0 allocations: 0 bytes)
# 1.912 μs (0 allocations: 0 bytes)
println("ISGate")
@btime benchmark($(state1), Clf.ISGate())
@btime benchmark($(state2), Clf.ISGate())

# 97.806 ns (0 allocations: 0 bytes)
# 1.995 μs (0 allocations: 0 bytes)
println("SXGate")
@btime benchmark($(state1), Clf.SGate())
@btime benchmark($(state2), Clf.SGate())

# 98.070 ns (0 allocations: 0 bytes)
# 1.925 μs (0 allocations: 0 bytes)
println("ISXGate")
@btime benchmark($(state1), Clf.ISGate())
@btime benchmark($(state2), Clf.ISGate())

# 3.151 μs (0 allocations: 0 bytes)
# 616.982 μs (0 allocations: 0 bytes)
println("CNOTGate")
@btime benchmark($(state1), Clf.CNOTGate())
@btime benchmark($(state2), Clf.CNOTGate())

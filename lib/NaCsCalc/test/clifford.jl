#!/usr/bin/julia

module TestClifford

using Test

using NaCsCalc

const Clf = NaCsCalc.Clifford

function test_gate(gate, inputs, outputs, sign)
    @test Clf.apply(gate, inputs..., false) === (outputs..., sign)
    @test Clf.apply(gate, inputs..., true) === (outputs..., !sign)
end

@testset "I" begin
    test_gate(Clf.IGate(), (false, false), (false, false), false)
    test_gate(Clf.IGate(), (true, false), (true, false), false)
    test_gate(Clf.IGate(), (true, true), (true, true), false)
    test_gate(Clf.IGate(), (false, true), (false, true), false)
end

@testset "X" begin
    test_gate(Clf.XGate(), (false, false), (false, false), false)
    test_gate(Clf.XGate(), (true, false), (true, false), false)
    test_gate(Clf.XGate(), (true, true), (true, true), true)
    test_gate(Clf.XGate(), (false, true), (false, true), true)

    test_gate(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
              (false, false), (false, false), false)
    test_gate(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
              (true, false), (true, false), false)
    test_gate(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
              (true, true), (true, true), true)
    test_gate(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
              (false, true), (false, true), true)
end

@testset "Y" begin
    test_gate(Clf.YGate(), (false, false), (false, false), false)
    test_gate(Clf.YGate(), (true, false), (true, false), true)
    test_gate(Clf.YGate(), (true, true), (true, true), false)
    test_gate(Clf.YGate(), (false, true), (false, true), true)
end

@testset "Z" begin
    test_gate(Clf.ZGate(), (false, false), (false, false), false)
    test_gate(Clf.ZGate(), (true, false), (true, false), true)
    test_gate(Clf.ZGate(), (true, true), (true, true), true)
    test_gate(Clf.ZGate(), (false, true), (false, true), false)

    test_gate(Clf.HGate() * Clf.XGate() * Clf.HGate(),
              (false, false), (false, false), false)
    test_gate(Clf.HGate() * Clf.XGate() * Clf.HGate(),
              (true, false), (true, false), true)
    test_gate(Clf.HGate() * Clf.XGate() * Clf.HGate(),
              (true, true), (true, true), true)
    test_gate(Clf.HGate() * Clf.XGate() * Clf.HGate(),
              (false, true), (false, true), false)
end

@testset "H" begin
    test_gate(Clf.HGate(), (false, false), (false, false), false)
    test_gate(Clf.HGate(), (true, false), (false, true), false)
    test_gate(Clf.HGate(), (true, true), (true, true), true)
    test_gate(Clf.HGate(), (false, true), (true, false), false)
end

@testset "S" begin
    test_gate(Clf.SGate(), (false, false), (false, false), false)
    test_gate(Clf.SGate(), (true, false), (true, true), false)
    test_gate(Clf.SGate(), (true, true), (true, false), true)
    test_gate(Clf.SGate(), (false, true), (false, true), false)

    test_gate(Clf.ISGate(), (false, false), (false, false), false)
    test_gate(Clf.ISGate(), (true, false), (true, true), true)
    test_gate(Clf.ISGate(), (true, true), (true, false), false)
    test_gate(Clf.ISGate(), (false, true), (false, true), false)
end

@testset "SX" begin
    test_gate(Clf.SXGate(), (false, false), (false, false), false)
    test_gate(Clf.SXGate(), (true, false), (true, false), false)
    test_gate(Clf.SXGate(), (true, true), (false, true), false)
    test_gate(Clf.SXGate(), (false, true), (true, true), true)

    test_gate(Clf.ISXGate(), (false, false), (false, false), false)
    test_gate(Clf.ISXGate(), (true, false), (true, false), false)
    test_gate(Clf.ISXGate(), (true, true), (false, true), true)
    test_gate(Clf.ISXGate(), (false, true), (true, true), false)
end

@testset "SY" begin
    test_gate(Clf.SYGate(), (false, false), (false, false), false)
    test_gate(Clf.SYGate(), (true, false), (false, true), true)
    test_gate(Clf.SYGate(), (true, true), (true, true), false)
    test_gate(Clf.SYGate(), (false, true), (true, false), false)

    test_gate(Clf.ISYGate(), (false, false), (false, false), false)
    test_gate(Clf.ISYGate(), (true, false), (false, true), false)
    test_gate(Clf.ISYGate(), (true, true), (true, true), false)
    test_gate(Clf.ISYGate(), (false, true), (true, false), true)
end

@testset "CNOT" begin
    test_gate(Clf.CNOTGate(), (false, false, false, false),
              (false, false, false, false), false)
    test_gate(Clf.CNOTGate(), (false, false, false, true),
              (false, true, false, true), false)
    test_gate(Clf.CNOTGate(), (false, false, true, false),
              (false, false, true, false), false)
    test_gate(Clf.CNOTGate(), (false, false, true, true),
              (false, true, true, true), false)
    test_gate(Clf.CNOTGate(), (false, true, false, false),
              (false, true, false, false), false)
    test_gate(Clf.CNOTGate(), (false, true, false, true),
              (false, false, false, true), false)
    test_gate(Clf.CNOTGate(), (false, true, true, false),
              (false, true, true, false), false)
    test_gate(Clf.CNOTGate(), (false, true, true, true),
              (false, false, true, true), false)

    test_gate(Clf.CNOTGate(), (true, false, false, false),
              (true, false, true, false), false)
    test_gate(Clf.CNOTGate(), (true, false, false, true),
              (true, true, true, true), true)
    test_gate(Clf.CNOTGate(), (true, false, true, false),
              (true, false, false, false), false)
    test_gate(Clf.CNOTGate(), (true, false, true, true),
              (true, true, false, true), false)
    test_gate(Clf.CNOTGate(), (true, true, false, false),
              (true, true, true, false), false)
    test_gate(Clf.CNOTGate(), (true, true, false, true),
              (true, false, true, true), false)
    test_gate(Clf.CNOTGate(), (true, true, true, false),
              (true, true, false, false), false)
    test_gate(Clf.CNOTGate(), (true, true, true, true),
              (true, false, false, true), true)
end

const inputs_1q = (0b1111_0000, 0b1100_1100, 0b1010_1010)
const gates_1q = [Clf.Clifford1Q{XZ[1],XZ[2],R[1],XZ[3],XZ[4],R[2]}()
                  for XZ in ((true, false, true, true),
                             (true, false, false, true),
                             (true, true, true, false),
                             (true, true, false, true),
                             (false, true, true, false),
                             (false, true, true, true)),
                      R in Iterators.product((false, true), (false, true))]

@testset "inv 1q" begin
    for gate in gates_1q
        @test gate * inv(gate) === Clf.IGate()
        @test Clf.apply(inv(gate), Clf.apply(gate, inputs_1q...)...) === inputs_1q
        @test Clf.apply(gate, Clf.apply(inv(gate), inputs_1q...)...) === inputs_1q
    end
end

@testset "prod 1q" begin
    for gate1 in gates_1q
        for gate2 in gates_1q
            @test(Clf.apply(gate2, Clf.apply(gate1, inputs_1q...)...) ===
                Clf.apply(gate1 * gate2, inputs_1q...))
            @test(Clf.apply(inv(gate1), Clf.apply(inv(gate2), inputs_1q...)...) ===
                Clf.apply(inv(gate1 * gate2), inputs_1q...))
            @test(Clf.apply(inv(gate2) * inv(gate1), inputs_1q...) ===
                Clf.apply(inv(gate1 * gate2), inputs_1q...))
        end
    end
end

function test_flip_base(::Type{SST}, n) where SST
    state = SST(n)

    res = Bool[]
    res0 = false
    for i in 1:n - 1
        v, det = Clf.measure_xs!(state, (i, i + 1))
        @test !det
        push!(res, v)
        res0 ⊻= v
    end
    # At this point, the stabilizers are all pairs of XX and the global Z
    v, det = Clf.measure_xs!(state, (1, n))
    @test det
    @test v == res0
    push!(res, v)
    for i in 1:n
        i2 = i == n ? 1 : i + 1
        v, det = Clf.measure_xs!(state, (i, i2))
        @test det
        @test v == res[i]
    end
    for i in 1:n
        for i2 in 1:n - 1
            if i2 >= i
                i2 += 1
            end
            v, det = Clf.measure_xs!(state, (i, i2))
            @test det
        end
    end
    v, det = Clf.measure_zs!(state, 1:n)
    @test det
    @test !v
    res_xn2, det = Clf.measure_xs!(state, 1:(n & ~1))
    @test det

    res = Bool[]
    res0 = false
    res1 = false
    for i in 1:n - 2
        v, det = Clf.measure_zs!(state, (i, i + 1))
        @test !det
        push!(res, v)
        res0 ⊻= v
        if isodd(i)
            res1 ⊻= v
        end
    end
    v, det = Clf.measure_zs!(state, 1:n)
    @test det
    @test !v
    # At this point, the first n-2 pairs of ZZ's are the stabilizer
    # and so is the global Z.
    # For even n, the combination of the ZZ's means that
    # Z(1)...Z(n-2) is also a stabilizer and therefore Z(n-1)Z(n) is as well.
    # The last stablizer in this case is X(1) ... X(n)
    # Otherwise, Z(n) is an stablizer
    # The last stablizer in this case is X(1) ... X(n - 1)
    if iseven(n)
        v, det = Clf.measure_zs!(state, (n - 1, n))
        @test det
        @test v == res1
        push!(res, v)
        res0 ⊻= v
        v, det = Clf.measure_xs!(state, 1:n)
        @test det
        @test v == res_xn2
        v, det = Clf.measure_zs!(state, (1, n))
        @test det
        @test v == res0
        push!(res, v)
    else
        v, det = Clf.measure_z!(state, n)
        @test det
        @test v == res1
        v, det = Clf.measure_xs!(state, 1:n - 1)
        @test det
        @test v == res_xn2

        v, det = Clf.measure_zs!(state, (n - 1, n))
        @test !det
        res0 ⊻= v
        push!(res, v)

        v, det = Clf.measure_zs!(state, (1, n))
        @test det
        @test v == res0
        push!(res, v)
    end
    for i in 1:n
        i2 = i == n ? 1 : i + 1
        v, det = Clf.measure_zs!(state, (i, i2))
        @test det
        @test v == res[i]
    end
    # At this point, all pairs of ZZ are the stabilizers
    for i in 1:n
        for i2 in 1:n - 1
            if i2 >= i
                i2 += 1
            end
            v, det = Clf.measure_zs!(state, (i, i2))
            @test det
        end
    end
    v, det = Clf.measure_zs!(state, 1:n)
    @test det
    @test !v
    if iseven(n)
        v, det = Clf.measure_xs!(state, 1:n)
        @test det
        @test v == res_xn2
    else
        for i in 1:n
            v, det = Clf.measure_z!(state, n)
            @test det
        end
    end
end

function test_flip_base_x(::Type{SST}, n) where SST
    state = SST(n)
    Clf.init_state_x!(state)

    res = Bool[]
    res0 = false
    for i in 1:n - 1
        v, det = Clf.measure_zs!(state, (i, i + 1))
        @test !det
        push!(res, v)
        res0 ⊻= v
    end
    # At this point, the stabilizers are all pairs of ZZ and the global X
    v, det = Clf.measure_zs!(state, (1, n))
    @test det
    @test v == res0
    push!(res, v)
    for i in 1:n
        i2 = i == n ? 1 : i + 1
        v, det = Clf.measure_zs!(state, (i, i2))
        @test det
        @test v == res[i]
    end
    for i in 1:n
        for i2 in 1:n - 1
            if i2 >= i
                i2 += 1
            end
            v, det = Clf.measure_zs!(state, (i, i2))
            @test det
        end
    end
    v, det = Clf.measure_xs!(state, 1:n)
    @test det
    @test !v
    res_zn2, det = Clf.measure_zs!(state, 1:(n & ~1))
    @test det

    res = Bool[]
    res0 = false
    res1 = false
    for i in 1:n - 2
        v, det = Clf.measure_xs!(state, (i, i + 1))
        @test !det
        push!(res, v)
        res0 ⊻= v
        if isodd(i)
            res1 ⊻= v
        end
    end
    v, det = Clf.measure_xs!(state, 1:n)
    @test det
    @test !v
    # At this point, the first n-2 pairs of XX's are the stabilizer
    # and so is the global X.
    # For even n, the combination of the XX's means that
    # X(1)...X(n-2) is also a stabilizer and therefore X(n-1)X(n) is as well.
    # The last stablizer in this case is Z(1) ... Z(n)
    # Otherwise, X(n) is an stablizer
    # The last stablizer in this case is Z(1) ... Z(n - 1)
    if iseven(n)
        v, det = Clf.measure_xs!(state, (n - 1, n))
        @test det
        @test v == res1
        push!(res, v)
        res0 ⊻= v
        v, det = Clf.measure_zs!(state, 1:n)
        @test det
        @test v == res_zn2
        v, det = Clf.measure_xs!(state, (1, n))
        @test det
        @test v == res0
        push!(res, v)
    else
        v, det = Clf.measure_x!(state, n)
        @test det
        @test v == res1
        v, det = Clf.measure_zs!(state, 1:n - 1)
        @test det
        @test v == res_zn2

        v, det = Clf.measure_xs!(state, (n - 1, n))
        @test !det
        res0 ⊻= v
        push!(res, v)

        v, det = Clf.measure_xs!(state, (1, n))
        @test det
        @test v == res0
        push!(res, v)
    end
    for i in 1:n
        i2 = i == n ? 1 : i + 1
        v, det = Clf.measure_xs!(state, (i, i2))
        @test det
        @test v == res[i]
    end
    # At this point, all pairs of XX are the stabilizers
    for i in 1:n
        for i2 in 1:n - 1
            if i2 >= i
                i2 += 1
            end
            v, det = Clf.measure_xs!(state, (i, i2))
            @test det
        end
    end
    v, det = Clf.measure_xs!(state, 1:n)
    @test det
    @test !v
    if iseven(n)
        v, det = Clf.measure_zs!(state, 1:n)
        @test det
        @test v == res_zn2
    else
        for i in 1:n
            v, det = Clf.measure_x!(state, n)
            @test det
        end
    end
end

function test_gap(::Type{SST}, gap, nqubit) where SST
    state = SST(nqubit)

    Clf.apply!(state, Clf.HGate(), 1)
    Clf.apply!(state, Clf.CNOTGate(), 1, gap + 1)
    Clf.apply!(state, Clf.HGate(), gap + 1)
    Clf.apply!(state, Clf.CNOTGate(), 1, 2gap + 1)
    Clf.apply!(state, Clf.HGate(), 2 * gap + 1)
    Clf.apply!(state, Clf.HGate(), 1)

    @test Clf.measure_xs!(state, (gap + 1, 2 * gap + 1)) == (false, true)

    meas, det = Clf.measure_z!(state, 1, force=true)
    @test !det

    @test Clf.measure_xs!(state, (gap + 1, 2 * gap + 1)) == (false, true)

    meas, det = Clf.measure_z!(state, 1)
    @test det

    @test Clf.measure_xs!(state, (gap + 1, 2 * gap + 1)) == (false, true)
end

function test_gap2(::Type{SST}, pos1, pos2, nqubit) where SST
    state = SST(nqubit)

    Clf.apply!(state, Clf.CNOTGate(), pos1, pos2)
    Clf.apply!(state, Clf.HGate(), pos1)
    Clf.apply!(state, Clf.CNOTGate(), pos1, pos2)
    if SST === Clf.InvStabilizerState
        inv_stab = Clf.get_inv_stabilizer(state, pos2, true)
        @test inv_stab.xs == [(i == pos1 || i == pos2) for i in 1:nqubit]
        @test inv_stab.zs == [(i == pos1 || i == pos2) for i in 1:nqubit]
        @test inv_stab.rs[] == true
    end
    Clf.apply!(state, Clf.SXGate(), pos2)
    Clf.apply!(state, Clf.SGate(), pos1)
    Clf.apply!(state, Clf.HGate(), pos1)
    Clf.apply!(state, Clf.CNOTGate(), pos1, pos2)

    @test Clf.measure_z!(state, pos2) == (false, true)
end

function apply_random_clifford!(cb, nbit)
    id = rand(1:32)
    n1 = rand(1:nbit)
    gate1q = get(gates_1q, id, Clf.IGate())
    if gate1q !== Clf.IGate()
        cb(gate1q, n1)
    else
        # CNOT
        n2 = rand(1:(nbit - 1))
        if n2 >= n1
            n2 += 1
        end
        cb(Clf.CNOTGate(), n1, n2)
    end
end

function test_deterministic_measure(::Type{SST}, nbit, ngates, nrep) where SST
    stabx = Vector{Bool}(undef, nbit)
    stabz = Vector{Bool}(undef, nbit)
    state = SST(nbit)
    str = Clf.PauliString{UInt64}(nbit)
    for i in 1:nrep
        if rand(Bool)
            Clf.init_state_z!(state)
            str.zs .= rand.(UInt64)
            str.xs .= zero(UInt64)
        else
            Clf.init_state_x!(state)
            str.zs .= zero(UInt64)
            str.xs .= rand.(UInt64)
        end
        str.rs[] = zero(UInt64)
        for _ in 1:ngates
            for i in 1:100
                apply_random_clifford!(nbit) do (args...)
                    Clf.apply!(state, args...)
                    Clf.apply!(str, args...)
                end
            end
            for j in 0:(Clf.nbits(typeof(str)) - 1)
                stabx .= ((str.xs .>> j) .& 1) .!= 0
                stabz .= ((str.zs .>> j) .& 1) .!= 0
                r = (str.rs[] >> j) & 1 != 0
                v, det = Clf.measure_paulis!(state, stabx, stabz)
                @test det
                @test r == v
            end
        end
    end
end

function rand_commuting_pauli!(strs, nbit)
    for str in strs
        str.zs .= rand.(UInt64)
        str.xs .= zero(UInt64)
        str.rs[] = zero(UInt64)
    end
    for i in 1:nbit * 20
        apply_random_clifford!(nbit) do (args...)
            for str in strs
                Clf.apply!(str, args...)
            end
        end
    end
end

function test_random_measure(::Type{SST}, nbit, ngates, nrep) where SST
    stabx = Vector{Bool}(undef, nbit)
    stabz = Vector{Bool}(undef, nbit)
    state = SST(nbit)
    strs = [Clf.PauliString{UInt64}(nbit) for i in 1:ngates]
    for i in 1:nrep
        if rand(Bool)
            Clf.init_state_z!(state)
        else
            Clf.init_state_x!(state)
        end
        rand_commuting_pauli!(strs, nbit)
        for j in 1:ngates
            for i in 1:20
                apply_random_clifford!(nbit) do (args...)
                    Clf.apply!(state, args...)
                    for str in strs
                        Clf.apply!(str, args...)
                    end
                end
            end
            for k in 1:(j - 1)
                str = strs[k]
                for b in 0:(Clf.nbits(typeof(str)) - 1)
                    stabx .= ((str.xs .>> b) .& 1) .!= 0
                    stabz .= ((str.zs .>> b) .& 1) .!= 0
                    r = (str.rs[] >> b) & 1 != 0
                    v, det = Clf.measure_paulis!(state, stabx, stabz)
                    @test det
                    @test r == v
                end
            end
            str = strs[j]
            for b in 0:(Clf.nbits(typeof(str)) - 1)
                stabx .= ((str.xs .>> b) .& 1) .!= 0
                stabz .= ((str.zs .>> b) .& 1) .!= 0
                r = (str.rs[] >> b) & 1 != 0
                v, det = Clf.measure_paulis!(state, stabx, stabz)
                str.rs[] = (v ? (str.rs[] | (UInt64(1) << b)) :
                    (str.rs[] & ~(UInt64(1) << b)))
            end
        end
    end
end

for (ss_name, SS_T) in [("StabilizerState", Clf.StabilizerState),
                        ("InvStabilizerState", Clf.InvStabilizerState)]
    @testset "Init ($ss_name)" begin
        for nbit in 1:3:2049
            state = SS_T(nbit)
            for i in 1:nbit
                @test Clf.measure_z!(state, i) == (false, true)
            end
            for v in (false, true)
                Clf.init_state_x!(state, v)
                for i in 1:nbit
                    @test Clf.measure_x!(state, i) == (v, true)
                end
                Clf.init_state_y!(state, v)
                for i in 1:nbit
                    @test Clf.measure_y!(state, i) == (v, true)
                end
                Clf.init_state_z!(state, v)
                for i in 1:nbit
                    @test Clf.measure_z!(state, i) == (v, true)
                end
            end
        end
    end

    @testset "stabilizer 1q ($ss_name)" begin
        state = SS_T(1)
        @test Clf.measure_z!(state, 1) === (false, true)
        Clf.apply!(state, Clf.XGate(), 1)
        @test Clf.measure_z!(state, 1) === (true, true)
        Clf.apply!(state, Clf.YGate(), 1)
        @test Clf.measure_z!(state, 1) === (false, true)
        Clf.apply!(state, Clf.ZGate(), 1)
        @test Clf.measure_z!(state, 1) === (false, true)

        v, det = Clf.measure_x!(state, 1)
        @test !det
        @test Clf.measure_x!(state, 1) === (v, true)
        @test Clf.measure_x!(state, 1) === (v, true)
        Clf.apply!(state, Clf.XGate(), 1)
        @test Clf.measure_x!(state, 1) === (v, true)
        Clf.apply!(state, Clf.YGate(), 1)
        @test Clf.measure_x!(state, 1) === (!v, true)
        Clf.apply!(state, Clf.ZGate(), 1)
        @test Clf.measure_x!(state, 1) === (v, true)

        v, det = Clf.measure_y!(state, 1)
        @test !det
        @test Clf.measure_y!(state, 1) === (v, true)
        @test Clf.measure_y!(state, 1) === (v, true)
        Clf.apply!(state, Clf.XGate(), 1)
        @test Clf.measure_y!(state, 1) === (!v, true)
        Clf.apply!(state, Clf.YGate(), 1)
        @test Clf.measure_y!(state, 1) === (!v, true)
        Clf.apply!(state, Clf.ZGate(), 1)
        @test Clf.measure_y!(state, 1) === (v, true)

        v, det = Clf.measure_z!(state, 1)
        @test !det
        @test Clf.measure_z!(state, 1) === (v, true)
        @test Clf.measure_z!(state, 1) === (v, true)
        if v
            Clf.apply!(state, Clf.XGate(), 1)
        end
        @test Clf.measure_z!(state, 1) === (false, true)

        @test Clf.measure_zs!(state, [1]) === (false, true)

        v, det = Clf.measure_xs!(state, [1])
        @test !det
        @test Clf.measure_x!(state, 1) === (v, true)
        @test Clf.measure_xs!(state, [1]) === (v, true)

        v, det = Clf.measure_ys!(state, [1])
        @test !det
        @test Clf.measure_y!(state, 1) === (v, true)
        @test Clf.measure_ys!(state, [1]) === (v, true)

        v, det = Clf.measure_zs!(state, [1])
        @test !det
        @test Clf.measure_z!(state, 1) === (v, true)
        @test Clf.measure_zs!(state, [1]) === (v, true)

        @test Clf.measure_paulis!(state, [false], [true]) === (v, true)

        v, det = Clf.measure_paulis!(state, [true], [false])
        @test !det
        @test Clf.measure_x!(state, 1) === (v, true)
        @test Clf.measure_paulis!(state, [true], [false]) === (v, true)

        v, det = Clf.measure_paulis!(state, [true], [true])
        @test !det
        @test Clf.measure_y!(state, 1) === (v, true)
        @test Clf.measure_paulis!(state, [true], [true]) === (v, true)

        v, det = Clf.measure_paulis!(state, [false], [true])
        @test !det
        @test Clf.measure_z!(state, 1) === (v, true)
        @test Clf.measure_paulis!(state, [false], [true]) === (v, true)
    end

    @testset "stabilizer 2q ($ss_name)" begin
        state = SS_T(2)
        @test Clf.measure_z!(state, 1) === (false, true)
        Clf.apply!(state, Clf.XGate(), 1)
        @test Clf.measure_z!(state, 1) === (true, true)
        Clf.apply!(state, Clf.YGate(), 1)
        @test Clf.measure_z!(state, 1) === (false, true)
        Clf.apply!(state, Clf.ZGate(), 1)
        @test Clf.measure_z!(state, 1) === (false, true)

        @test Clf.measure_z!(state, 2) === (false, true)
        Clf.apply!(state, Clf.XGate(), 2)
        @test Clf.measure_z!(state, 2) === (true, true)
        Clf.apply!(state, Clf.YGate(), 2)
        @test Clf.measure_z!(state, 2) === (false, true)
        Clf.apply!(state, Clf.ZGate(), 2)
        @test Clf.measure_z!(state, 2) === (false, true)

        v, det = Clf.measure_x!(state, 1)
        @test !det
        @test Clf.measure_x!(state, 1) === (v, true)
        @test Clf.measure_x!(state, 1) === (v, true)
        Clf.apply!(state, Clf.XGate(), 1)
        @test Clf.measure_x!(state, 1) === (v, true)
        Clf.apply!(state, Clf.YGate(), 1)
        @test Clf.measure_x!(state, 1) === (!v, true)
        Clf.apply!(state, Clf.ZGate(), 1)
        @test Clf.measure_x!(state, 1) === (v, true)

        v, det = Clf.measure_y!(state, 1)
        @test !det
        @test Clf.measure_y!(state, 1) === (v, true)
        @test Clf.measure_y!(state, 1) === (v, true)
        Clf.apply!(state, Clf.XGate(), 1)
        @test Clf.measure_y!(state, 1) === (!v, true)
        Clf.apply!(state, Clf.YGate(), 1)
        @test Clf.measure_y!(state, 1) === (!v, true)
        Clf.apply!(state, Clf.ZGate(), 1)
        @test Clf.measure_y!(state, 1) === (v, true)

        v, det = Clf.measure_z!(state, 1)
        @test !det
        @test Clf.measure_z!(state, 1) === (v, true)
        @test Clf.measure_z!(state, 1) === (v, true)
        if v
            Clf.apply!(state, Clf.XGate(), 1)
        end
        @test Clf.measure_z!(state, 1) === (false, true)

        v, det = Clf.measure_x!(state, 2)
        @test !det
        @test Clf.measure_x!(state, 2) === (v, true)
        @test Clf.measure_x!(state, 2) === (v, true)
        Clf.apply!(state, Clf.XGate(), 2)
        @test Clf.measure_x!(state, 2) === (v, true)
        Clf.apply!(state, Clf.YGate(), 2)
        @test Clf.measure_x!(state, 2) === (!v, true)
        Clf.apply!(state, Clf.ZGate(), 2)
        @test Clf.measure_x!(state, 2) === (v, true)

        v, det = Clf.measure_y!(state, 2)
        @test !det
        @test Clf.measure_y!(state, 2) === (v, true)
        @test Clf.measure_y!(state, 2) === (v, true)
        Clf.apply!(state, Clf.XGate(), 2)
        @test Clf.measure_y!(state, 2) === (!v, true)
        Clf.apply!(state, Clf.YGate(), 2)
        @test Clf.measure_y!(state, 2) === (!v, true)
        Clf.apply!(state, Clf.ZGate(), 2)
        @test Clf.measure_y!(state, 2) === (v, true)

        v, det = Clf.measure_z!(state, 2)
        @test !det
        @test Clf.measure_z!(state, 2) === (v, true)
        @test Clf.measure_z!(state, 2) === (v, true)
        if v
            Clf.apply!(state, Clf.XGate(), 2)
        end
        @test Clf.measure_z!(state, 2) === (false, true)

        Clf.apply!(state, Clf.HGate(), 1)
        v, det = Clf.measure_zs!(state, [1, 2])
        @test !det
        @test Clf.measure_z!(state, 1) === (v, true)
        @test Clf.measure_z!(state, 2) === (false, true)

        Clf.apply!(state, Clf.HGate(), 1)
        v2, det = Clf.measure_xs!(state, [1, 2])
        @test !det
        @test Clf.measure_x!(state, 1) === (v, true)
        @test Clf.measure_x!(state, 2) === (v2 ⊻ v, true)

        v, det = Clf.measure_ys!(state, [1, 2])
        @test !det
        v2, det = Clf.measure_y!(state, 1)
        @test Clf.measure_y!(state, 1) === (v2, true)
        @test Clf.measure_y!(state, 2) === (v2 ⊻ v, true)

        Clf.init_state_z!(state)
        @test Clf.measure_z!(state, 1) === (false, true)
        @test Clf.measure_z!(state, 2) === (false, true)

        Clf.apply!(state, Clf.HGate(), 1)
        Clf.apply!(state, Clf.HGate(), 2)
        @test Clf.measure_x!(state, 1) === (false, true)
        @test Clf.measure_x!(state, 2) === (false, true)

        v, det = Clf.measure_zs!(state, [1, 2])
        @test !det
        v2, det = Clf.measure_z!(state, 1)
        @test !det
        @test Clf.measure_z!(state, 1) === (v2, true)
        @test Clf.measure_z!(state, 2) === (v2 ⊻ v, true)

        v, det = Clf.measure_paulis!(state, [true, true], [false, true])
        @test !det
        v2, det = Clf.measure_x!(state, 1)
        @test !det
        @test Clf.measure_x!(state, 1) === (v2, true)
        @test Clf.measure_y!(state, 2) === (v2 ⊻ v, true)
    end

    @testset "flip base ($ss_name)" begin
        test_flip_base(SS_T, 4)
        test_flip_base(SS_T, 5)
        test_flip_base(SS_T, 8)
        test_flip_base(SS_T, 16)
        test_flip_base(SS_T, 17)
        test_flip_base(SS_T, 31)
        test_flip_base(SS_T, 32)
        test_flip_base(SS_T, 64)
        test_flip_base(SS_T, 65)
        test_flip_base(SS_T, 70)
        test_flip_base(SS_T, 80)
        test_flip_base(SS_T, 90)
        test_flip_base(SS_T, 100)
        test_flip_base(SS_T, 110)
        test_flip_base(SS_T, 120)
        test_flip_base(SS_T, 127)
        test_flip_base(SS_T, 128)
        test_flip_base(SS_T, 129)

        test_flip_base_x(SS_T, 4)
        test_flip_base_x(SS_T, 5)
        test_flip_base_x(SS_T, 16)
        test_flip_base_x(SS_T, 17)
        test_flip_base_x(SS_T, 31)
        test_flip_base_x(SS_T, 32)
        test_flip_base_x(SS_T, 64)
        test_flip_base_x(SS_T, 65)
        test_flip_base_x(SS_T, 70)
        test_flip_base_x(SS_T, 80)
        test_flip_base_x(SS_T, 90)
        test_flip_base_x(SS_T, 100)
        test_flip_base_x(SS_T, 110)
        test_flip_base_x(SS_T, 120)
        test_flip_base_x(SS_T, 127)
        test_flip_base_x(SS_T, 128)
        test_flip_base_x(SS_T, 129)
    end

    @testset "measure sign ($ss_name)" begin
        state = SS_T(3)
        @test Clf.measure_z!(state, 3) === (false, true)
        Clf.apply!(state, Clf.CNOTGate(), 2, 1)
        Clf.apply!(state, Clf.CNOTGate(), 1, 3)
        Clf.apply!(state, Clf.HGate(), 1)
        Clf.apply!(state, Clf.CNOTGate(), 1, 2)
        @test Clf.measure_z!(state, 3) === (false, true)
    end

    @testset "gap ($ss_name)" begin
        for gap in 1:128
            for nqubit in [3:140; 240:260; 500:530; 1000:1050; 2020:2070]
                if nqubit < gap * 2 + 1
                    continue
                end
                test_gap(SS_T, gap, nqubit)
            end
        end
    end

    @testset "gap2 ($ss_name)" begin
        for nqubit in [5; 63:65; 127:129; 250:2:262;
                       504:4:524; 1008:8:1040; 2016:16:2080]
            step = nqubit >= 1000 ? 8 : (nqubit > 500 ? 4 : 2)
            for pos1 in 1:step:nqubit
                for pos2 in 1:step:(nqubit - 1)
                    if pos2 >= pos1
                        pos2 += 1
                    end
                    test_gap2(SS_T, pos1, pos2, nqubit)
                end
            end
        end
    end

    @testset "random measurement ($ss_name)" begin
        test_random_measure(SS_T, 10, 80, 100)
        test_random_measure(SS_T, 31, 80, 30)
        test_random_measure(SS_T, 32, 80, 30)
        test_random_measure(SS_T, 33, 80, 30)
        test_random_measure(SS_T, 63, 80, 30)
        test_random_measure(SS_T, 64, 80, 30)
        test_random_measure(SS_T, 65, 80, 30)
        test_random_measure(SS_T, 97, 80, 25)
        test_random_measure(SS_T, 128, 80, 20)
    end

    @testset "deterministic measurement ($ss_name)" begin
        test_deterministic_measure(SS_T, 10, 2000, 100)
        test_deterministic_measure(SS_T, 31, 2000, 30)
        test_deterministic_measure(SS_T, 32, 2000, 30)
        test_deterministic_measure(SS_T, 33, 2000, 30)
        test_deterministic_measure(SS_T, 63, 2000, 30)
        test_deterministic_measure(SS_T, 64, 2000, 30)
        test_deterministic_measure(SS_T, 65, 2000, 30)
        test_deterministic_measure(SS_T, 97, 2000, 25)
        test_deterministic_measure(SS_T, 128, 2000, 20)
    end
end

end

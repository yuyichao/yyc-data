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

const inputs_1q = Iterators.product((false, true), (false, true), (false, true))

@testset "inv" begin
    for gate in [Clf.IGate(), Clf.XGate(), Clf.YGate(), Clf.ZGate(), Clf.HGate(),
                 Clf.SGate(), Clf.ISGate(), Clf.SXGate(), Clf.ISXGate()]
        for input in inputs_1q
            @test Clf.apply(inv(gate), Clf.apply(gate, input...)...) === input
            @test Clf.apply(gate, Clf.apply(inv(gate), input...)...) === input
        end
    end
end

@testset "prod" begin
    for gate1 in [Clf.IGate(), Clf.XGate(), Clf.YGate(), Clf.ZGate(), Clf.HGate(),
                  Clf.SGate(), Clf.ISGate(), Clf.SXGate(), Clf.ISXGate()]
        for gate2 in [Clf.IGate(), Clf.XGate(), Clf.YGate(), Clf.ZGate(), Clf.HGate(),
                      Clf.SGate(), Clf.ISGate(), Clf.SXGate(), Clf.ISXGate()]
            for input in inputs_1q
                @test(Clf.apply(gate2, Clf.apply(gate1, input...)...) ===
                    Clf.apply(gate1 * gate2, input...))
                @test(Clf.apply(inv(gate1), Clf.apply(inv(gate2), input...)...) ===
                    Clf.apply(inv(gate1 * gate2), input...))
                @test(Clf.apply(inv(gate2) * inv(gate1), input...) ===
                    Clf.apply(inv(gate1 * gate2), input...))
            end
        end
    end
end

@testset "stablizer 1q" begin
    state = Clf.StabilizerState(1)
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

@testset "stablizer 2q" begin
    state = Clf.StabilizerState(2)
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

function test_flip_base(n)
    state = Clf.StabilizerState(n)

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

function test_flip_base_x(n)
    state = Clf.StabilizerState(n)
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

@testset "flip base" begin
    test_flip_base(4)
    test_flip_base(5)
    test_flip_base(16)
    test_flip_base(17)
    test_flip_base(31)
    test_flip_base(32)
    test_flip_base(64)
    test_flip_base(65)
    test_flip_base(70)
    test_flip_base(80)
    test_flip_base(90)
    test_flip_base(100)
    test_flip_base(110)
    test_flip_base(120)
    test_flip_base(127)
    test_flip_base(128)
    test_flip_base(129)

    test_flip_base_x(4)
    test_flip_base_x(5)
    test_flip_base_x(16)
    test_flip_base_x(17)
    test_flip_base_x(31)
    test_flip_base_x(32)
    test_flip_base_x(64)
    test_flip_base_x(65)
    test_flip_base_x(70)
    test_flip_base_x(80)
    test_flip_base_x(90)
    test_flip_base_x(100)
    test_flip_base_x(110)
    test_flip_base_x(120)
    test_flip_base_x(127)
    test_flip_base_x(128)
    test_flip_base_x(129)
end

@testset "measure sign" begin
    state = Clf.StabilizerState(3)
    @test Clf.measure_z!(state, 3) === (false, true)
    Clf.apply!(state, Clf.CNOTGate(), 2, 1)
    Clf.apply!(state, Clf.CNOTGate(), 1, 3)
    Clf.apply!(state, Clf.HGate(), 1)
    Clf.apply!(state, Clf.CNOTGate(), 1, 2)
    @test Clf.measure_z!(state, 3) === (false, true)
end

end

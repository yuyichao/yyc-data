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

end

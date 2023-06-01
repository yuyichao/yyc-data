#!/usr/bin/julia

module TestClifford

using Test

using NaCsCalc

const Clf = NaCsCalc.Clifford

function test_gate_1q(gate, inputs, outputs, sign)
    @test Clf.apply(gate, inputs..., false) === (outputs..., sign)
    @test Clf.apply(gate, inputs..., true) === (outputs..., !sign)
end

@testset "X" begin
    test_gate_1q(Clf.XGate(), (false, false), (false, false), false)
    test_gate_1q(Clf.XGate(), (true, false), (true, false), false)
    test_gate_1q(Clf.XGate(), (true, true), (true, true), true)
    test_gate_1q(Clf.XGate(), (false, true), (false, true), true)

    test_gate_1q(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
                 (false, false), (false, false), false)
    test_gate_1q(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
                 (true, false), (true, false), false)
    test_gate_1q(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
                 (true, true), (true, true), true)
    test_gate_1q(Clf.HGate() * Clf.ZGate() * Clf.HGate(),
                 (false, true), (false, true), true)
end

@testset "Y" begin
    test_gate_1q(Clf.YGate(), (false, false), (false, false), false)
    test_gate_1q(Clf.YGate(), (true, false), (true, false), true)
    test_gate_1q(Clf.YGate(), (true, true), (true, true), false)
    test_gate_1q(Clf.YGate(), (false, true), (false, true), true)
end

@testset "Z" begin
    test_gate_1q(Clf.ZGate(), (false, false), (false, false), false)
    test_gate_1q(Clf.ZGate(), (true, false), (true, false), true)
    test_gate_1q(Clf.ZGate(), (true, true), (true, true), true)
    test_gate_1q(Clf.ZGate(), (false, true), (false, true), false)

    test_gate_1q(Clf.HGate() * Clf.XGate() * Clf.HGate(),
                 (false, false), (false, false), false)
    test_gate_1q(Clf.HGate() * Clf.XGate() * Clf.HGate(),
                 (true, false), (true, false), true)
    test_gate_1q(Clf.HGate() * Clf.XGate() * Clf.HGate(),
                 (true, true), (true, true), true)
    test_gate_1q(Clf.HGate() * Clf.XGate() * Clf.HGate(),
                 (false, true), (false, true), false)
end

@testset "H" begin
    test_gate_1q(Clf.HGate(), (false, false), (false, false), false)
    test_gate_1q(Clf.HGate(), (true, false), (false, true), false)
    test_gate_1q(Clf.HGate(), (true, true), (true, true), true)
    test_gate_1q(Clf.HGate(), (false, true), (true, false), false)
end

@testset "S" begin
    test_gate_1q(Clf.SGate(), (false, false), (false, false), false)
    test_gate_1q(Clf.SGate(), (true, false), (true, true), false)
    test_gate_1q(Clf.SGate(), (true, true), (true, false), true)
    test_gate_1q(Clf.SGate(), (false, true), (false, true), false)
end

@testset "SX" begin
    test_gate_1q(Clf.SXGate(), (false, false), (false, false), false)
    test_gate_1q(Clf.SXGate(), (true, false), (true, false), false)
    test_gate_1q(Clf.SXGate(), (true, true), (false, true), false)
    test_gate_1q(Clf.SXGate(), (false, true), (true, true), true)
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

end

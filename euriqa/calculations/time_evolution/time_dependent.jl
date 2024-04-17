#!/usr/bin/julia

using QuantumOptics
using LinearAlgebra

struct LabBasisData{B,K,O} <: Function
    Ω0_2::Float64
    dΩ_2::Float64
    δ0::Float64
    dδ_2::Float64
    tlen::Float64

    basis::B
    ψ0::K
    op::O
    function LabBasisData(Ω0, Ω1, δ0, δ1, tlen)
        basis = SpinBasis(1//2)
        op_data = zeros(ComplexF64, 2, 2)
        op = Operator(basis, op_data)
        ψ0 = spindown(basis)
        Ω0_2 = Ω0 / 2
        dΩ_2 = (Ω1 - Ω0) / tlen / 2
        dδ_2 = (δ1 - δ0) / tlen / 2
        return new{typeof(basis),typeof(ψ0),typeof(op)}(Ω0_2, dΩ_2, δ0, dδ_2,
                                                         tlen, basis, ψ0, op)
    end
end

function (data::LabBasisData)(t, ψ)
    op_data = data.op.data
    # δ = δ0 + (δ1 - δ0) * t / tlen
    # φ = ∫δdt = δ0 * t + (δ1 - δ0) * t^2 / 2 / tlen
    φ = t * muladd(data.dδ_2, t, data.δ0)
    Ω_2 = muladd(data.dΩ_2, t, data.Ω0_2)
    M = Ω_2 * cis(φ)
    op_data[1, 2] = M
    op_data[2, 1] = conj(M)
    return data.op
end

struct LabBasisCalc{LD}
    data::LD
    function LabBasisCalc(args...)
        ld = LabBasisData(args...)
        return new{typeof(ld)}(ld)
    end
end

function evolve(calc::LabBasisCalc, npoints=1001)
    calc_data = calc.data
    ts, ψs = timeevolution.schroedinger_dynamic(range(0, calc_data.tlen, npoints),
                                                 calc_data.ψ0, calc_data)
    return ts, [ψ.data[i] for i in 1:2, ψ in ψs]
end

struct RotBasisData{B,K,O} <: Function
    Ω0_2::Float64
    dΩ_2::Float64
    δ0_2::Float64
    dδ_2::Float64
    tlen::Float64
    basis::B
    ψ0::K
    op::O
    function RotBasisData(Ω0, Ω1, δ0, δ1, tlen)
        basis = SpinBasis(1//2)
        op_data = zeros(ComplexF64, 2, 2)
        op = Operator(basis, op_data)
        ψ0 = spindown(basis)
        Ω0_2 = Ω0 / 2
        dΩ_2 = (Ω1 - Ω0) / tlen / 2
        dδ_2 = (δ1 - δ0) / tlen / 2
        return new{typeof(basis),typeof(ψ0),typeof(op)}(Ω0_2, dΩ_2, δ0 / 2, dδ_2, tlen,
                                                         basis, ψ0, op)
    end
end

function (data::RotBasisData)(t, ψ)
    op_data = data.op.data
    δ_2 = muladd(data.dδ_2, t, data.δ0_2)
    Ω_2 = muladd(data.dΩ_2, t, data.Ω0_2)
    op_data[1, 1] = δ_2
    op_data[2, 2] = -δ_2
    op_data[1, 2] = Ω_2
    op_data[2, 1] = Ω_2
    return data.op
end

struct RotBasisCalc{LD}
    data::LD
    function RotBasisCalc(args...)
        ld = RotBasisData(args...)
        return new{typeof(ld)}(ld)
    end
end

function evolve(calc::RotBasisCalc, npoints=1001)
    calc_data = calc.data
    ts, ψs = timeevolution.schroedinger_dynamic(range(0, calc_data.tlen, npoints),
                                                 calc_data.ψ0, calc_data)
    data = [ψ.data[i] for i in 1:2, ψ in ψs]
    δ0_2 = calc_data.δ0_2
    dδ_4 = calc_data.dδ_2 / 2
    tlen = calc_data.tlen
    for (i, t) in enumerate(ts)
        φ_2 = t * muladd(dδ_4, t, δ0_2)
        p = cis(φ_2)
        data[1, i] *= p
        data[2, i] *= conj(p)
    end
    return ts, data
end

struct MagnusCalc
    A0::Matrix{ComplexF64}
    A1::Matrix{ComplexF64}
    A01::Matrix{ComplexF64}
    A011::Matrix{ComplexF64}
    A0100::Matrix{ComplexF64}
    A0101::Matrix{ComplexF64}
    A0110::Matrix{ComplexF64}
    A0111::Matrix{ComplexF64}

    Abuff::Matrix{ComplexF64}
    ψbuff1::Vector{ComplexF64}
    ψbuff2::Vector{ComplexF64}
    ψ0::Vector{ComplexF64}

    δ0_2::Float64
    dδ_4::Float64
    tlen::Float64

    function MagnusCalc(Ω0, Ω1, δ0, δ1, tlen)
        H0 = [δ0 / 2 Ω0 / 2
              Ω0 / 2 -δ0 / 2]
        H1 = [(δ1 - δ0) / 2 (Ω1 - Ω0) / 2
              (Ω1 - Ω0) / 2 -(δ1 - δ0) / 2]
        H1 ./= tlen

        A0 = H0 .* -im
        A1 = H1 .* -im

        A01 = A0 * A1 - A1 * A0

        A010 = A01 * A0 - A0 * A01
        A011 = A01 * A1 - A1 * A01

        A0100 = A010 * A0 - A0 * A010
        A0101 = A010 * A1 - A1 * A010
        A0110 = A011 * A0 - A0 * A011
        A0111 = A011 * A1 - A1 * A011

        ψ0 = [0, 1]

        return new(A0, A1, A01, A011, A0100, A0101, A0110, A0111,
                   zeros(ComplexF64, 2, 2), zeros(ComplexF64, 2), zeros(ComplexF64, 2),
                   ψ0, δ0 / 2, (δ1 - δ0) / tlen / 4, tlen)
    end
end

function evolve(calc::MagnusCalc, npoints=1001)
    tlen = calc.tlen

    ts = [range(0, tlen, npoints);]
    data = Matrix{ComplexF64}(undef, 2, npoints)

    A0 = calc.A0
    A1 = calc.A1
    A01 = calc.A01
    A011 = calc.A011
    A0100 = calc.A0100
    A0101 = calc.A0101
    A0110 = calc.A0110
    A0111 = calc.A0111

    Abuff = calc.Abuff
    ψbuff1 = calc.ψbuff1
    ψbuff2 = calc.ψbuff2

    data[:, 1] .= calc.ψ0
    ψbuff1 .= calc.ψ0
    δ0_2 = calc.δ0_2
    dδ_4 = calc.dδ_4

    @inbounds for i in 2:npoints
        t1 = ts[i - 1]
        t2 = ts[i]
        dt = t2 - t1
        Abuff .= (dt .* (A0 .+ (t1 + t2) / 2 .* A1) .- dt^3 / 12 .* A01 .+
            dt^5 / 240 .* (A011 .+ 1 / 3 .* (
                A0100 .+ (t1 + 5 * t2) / 6 .* A0101 .+
                    (5 * t1 + t2) / 6 .* A0110 .+
                    (t1^2 + 5 * t1 * t2 + t2^2) / 7 .* A0111)))
        LinearAlgebra.mul!(ψbuff2, LinearAlgebra.exp!(Abuff), ψbuff1)

        φ_2 = t2 * muladd(dδ_4, t2, δ0_2)
        p = cis(φ_2)
        data[1, i] = ψbuff2[1] * p
        data[2, i] = ψbuff2[2] * conj(p)

        ψbuff1, ψbuff2 = ψbuff2, ψbuff1
    end
    return ts, data
end

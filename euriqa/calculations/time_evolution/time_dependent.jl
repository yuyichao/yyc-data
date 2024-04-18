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

@inline tcross(v1, v2) = (muladd(v1[2], v2[3], -v1[3] * v2[2]),
                          muladd(v1[3], v2[1], -v1[1] * v2[3]),
                          muladd(v1[1], v2[2], - v1[2] * v2[1]))
@inline function expiv_pauli(b, v)
    len2 = muladd(v[1], v[1], muladd(v[2], v[2], abs2(v[3])))
    if len2 == 0
        return b
    end
    len = sqrt(len2)
    v = v ./ len
    s, c = @inline sincos(len)
    @inbounds begin
        X = s * v[1]
        Y = s * v[2]
        Z = s * v[3]

        U11 = complex(c, Z)
        U22 = complex(c, -Z)
        U12 = complex(Y, X)
        U21 = complex(-Y, X)

        return (muladd(U11, b[1], U12 * b[2]), muladd(U21, b[1], U22 * b[2]))
    end
end

struct MagnusAnalyticCalc
    h0::NTuple{3,Float64}
    h1::NTuple{3,Float64}

    h0ch1::NTuple{3,Float64}
    h0dh0::Float64
    h0dh1::Float64
    h1dh1::Float64

    Abuff::Matrix{ComplexF64}
    ψbuff1::Vector{ComplexF64}
    ψbuff2::Vector{ComplexF64}
    ψ0::Vector{ComplexF64}

    δ0_2::Float64
    dδ_4::Float64
    tlen::Float64

    function MagnusAnalyticCalc(Ω0, Ω1, δ0, δ1, tlen)
        h0 = (Ω0 / 2, 0.0, δ0 / 2)
        h1 = ((Ω1 - Ω0) / 2, 0.0, (δ1 - δ0) / 2) ./ tlen
        h0ch1 = tcross(h0, h1)

        h0dh0 = dot(h0, h0)
        h0dh1 = dot(h0, h1)
        h1dh1 = dot(h1, h1)

        ψ0 = ComplexF64[0, 1]

        return new(h0, h1, h0ch1, h0dh0, h0dh1, h1dh1,
                   zeros(ComplexF64, 2, 2), zeros(ComplexF64, 2), zeros(ComplexF64, 2),
                   ψ0, δ0 / 2, (δ1 - δ0) / tlen / 4, tlen)
    end
end

function evolve(calc::MagnusAnalyticCalc, npoints=1001)
    tlen = calc.tlen

    ts = [range(0, tlen, npoints);]
    data = Matrix{ComplexF64}(undef, 2, npoints)

    h0 = calc.h0
    h1 = calc.h1
    h0ch1 = calc.h0ch1
    h0dh0 = calc.h0dh0
    h0dh1 = calc.h0dh1
    h1dh1 = calc.h1dh1

    Abuff = calc.Abuff
    ψbuff1 = calc.ψbuff1
    ψbuff2 = calc.ψbuff2

    @inbounds begin
        data[:, 1] .= calc.ψ0
        ψ = (calc.ψ0[1], calc.ψ0[2])
    end
    δ0_2 = calc.δ0_2
    dδ_4 = calc.dδ_4

    @inbounds for i in 2:npoints
        t1 = ts[i - 1]
        t2 = ts[i]
        dt = t2 - t1
        dt2 = dt^2
        dt4 = dt2^2

        A = t1 + t2
        B = muladd(t1, t1, muladd(5 * t1, t2, t2^2))
        C = muladd(B, h1dh1 / 7, muladd(A, h0dh1, h0dh0))

        h = -dt .* (muladd.(dt4 / 60, muladd.(h0dh1, h1, .-h1dh1 .* h0),
                            muladd.(A / 2, h1, h0))
                    .- muladd(dt4 / 90, C, dt2 / 6) .* h0ch1)
        ψ = expiv_pauli(ψ, h)

        φ_2 = t2 * muladd(dδ_4, t2, δ0_2)
        p = @inline cis(φ_2)
        data[1, i] = ψ[1] * p
        data[2, i] = ψ[2] * conj(p)
    end
    return ts, data
end

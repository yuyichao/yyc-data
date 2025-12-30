#!/usr/bin/julia

using Test

include("sparse.jl")

function test_param(f!, C, X, A, α, β)
    X_orig = copy(X)
    A_orig = copy(A)

    C1 = copy(C)
    mul!(C1, X, A, α, β)

    C2 = copy(C)
    f!(C2, X, A, α, β)

    @test C1 ≈ C2
end

function test_matrix(f!, C, X, A)
    v0 = zero(eltype(C))
    v1 = one(eltype(C))
    v05 = v1 / 2

    vs = Any[false, true, v0, v05, v1]
    for α in vs
        for β in vs
            test_param(f!, C, X, A, α, β)
        end
    end
end

function test_type(f!, ElType, sz)
    for _ in 1:100
        test_matrix(f!, rand(ElType, sz, sz), rand(ElType, sz, sz),
                    sprand(ElType, sz, sz, 0.25))
        test_matrix(f!, rand(ElType, sz, sz), rand(ElType, sz, sz),
                    sprand(ElType, sz, sz, 0.05))
        test_matrix(f!, rand(ElType, sz, sz), rand(ElType, sz, sz),
                    sprand(ElType, sz, sz, 0.01))
        test_matrix(f!, rand(ElType, sz, sz), rand(ElType, sz, sz),
                    spzeros(ElType, sz, sz))
        test_matrix(f!, rand(ElType, sz, sz), rand(ElType, sz, sz),
                    sparse(1:sz, 1:sz, rand(ElType, sz)))
    end
end

function test_f(f!)
    for sz in [1, 3, 10, 30, 100]
        for ElType in [Float32, Float64, Int32, Int64, ComplexF32, ComplexF64, BigFloat]
            test_type(f!, ElType, sz)
        end
    end
end

@testset "Testing $(f!)" for f! in (spmul_orig, spmul_view, spmul_muladd, spmul_split)
    test_f(f!)
end

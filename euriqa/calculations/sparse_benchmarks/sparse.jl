#!/usr/bin/julia

using LinearAlgebra
using SparseArrays
using SparseArrays: DenseMatrixUnion, SparseMatrixCSCUnion2

using BenchmarkTools

function spmul1(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    @inbounds for col in axes(A,2), k in nzrange(A, col)
        Aiα = nzv[k] * α
        rvk = rv[k]
        @simd for multivec_row in axes(X,1)
            C[multivec_row, col] += X[multivec_row, rvk] * Aiα
        end
    end
    C
end

function spmul2(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    @inbounds for col in axes(A,2)
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            @simd for multivec_row in axes(X,1)
                C_col[multivec_row] += X_col[multivec_row] * Aiα
            end
        end
    end
    C
end

function spmul3(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β_one = isone(β)
    β_zero = iszero(β)
    @inbounds for col in axes(A,2)
        filled = β_one
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            if filled
                @simd for multivec_row in axes(X,1)
                    C_col[multivec_row] += X_col[multivec_row] * Aiα
                end
            elseif β_zero
                @simd for multivec_row in axes(X,1)
                    C_col[multivec_row] = X_col[multivec_row] * Aiα
                end
            else
                @simd for multivec_row in axes(X,1)
                    C_col[multivec_row] = C_col[multivec_row] * β + X_col[multivec_row] * Aiα
                end
            end
            filled = true
        end
        if !filled
            if β_zero
                fill!(C_col, zero(eltype(C)))
            else
                rmul!(C_col, β)
            end
        end
    end
    C
end

function spmul4(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β_one = isone(β)
    β_zero = iszero(β)
    @inbounds for col in axes(A,2)
        filled = β_one
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            if filled
                @simd for multivec_row in axes(X,1)
                    C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
                                                 C_col[multivec_row])
                end
            elseif β_zero
                @simd for multivec_row in axes(X,1)
                    C_col[multivec_row] = X_col[multivec_row] * Aiα
                end
            else
                @simd for multivec_row in axes(X,1)
                    C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
                                                 C_col[multivec_row] * β)
                end
            end
            filled = true
        end
        if !filled
            if β_zero
                fill!(C_col, zero(eltype(C)))
            else
                rmul!(C_col, β)
            end
        end
    end
    C
end

function bench_sparse(C, X, A)
    println("  false")
    @btime spmul1($C, $X, $A, true, false)
    @btime spmul2($C, $X, $A, true, false)
    @btime spmul3($C, $X, $A, true, false)
    @btime spmul4($C, $X, $A, true, false)

    println("  true")
    @btime spmul1($C, $X, $A, true, true)
    @btime spmul2($C, $X, $A, true, true)
    @btime spmul3($C, $X, $A, true, true)
    @btime spmul4($C, $X, $A, true, true)

    return
end

function bench_sparse_type(ElType, sz)
    println(" 25%")
    bench_sparse(zeros(ElType, sz, sz), zeros(ElType, sz, sz),
                 sprand(ElType, sz, sz, 0.25))
    println(" 5%")
    bench_sparse(zeros(ElType, sz, sz), zeros(ElType, sz, sz),
                 sprand(ElType, sz, sz, 0.05))
    println(" 1%")
    bench_sparse(zeros(ElType, sz, sz), zeros(ElType, sz, sz),
                 sprand(ElType, sz, sz, 0.01))
    println(" diag")
    bench_sparse(zeros(ElType, sz, sz), zeros(ElType, sz, sz),
                 sparse(1:sz, 1:sz, ones(ElType, sz)))
end

println("ComplexF64 10x10")
bench_sparse_type(ComplexF64, 10)
println("ComplexF64 100x100")
bench_sparse_type(ComplexF64, 100)
println("ComplexF64 1000x1000")
bench_sparse_type(ComplexF64, 1000)

println("Float64 10x10")
bench_sparse_type(Float64, 10)
println("Float64 100x100")
bench_sparse_type(Float64, 100)
println("Float64 1000x1000")
bench_sparse_type(Float64, 1000)

println("BigFloat 10x10")
bench_sparse_type(BigFloat, 10)
println("BigFloat 100x100")
bench_sparse_type(BigFloat, 100)
# println("BigFloat 1000x1000")
# bench_sparse_type(BigFloat, 1000)

println("Int64 10x10")
bench_sparse_type(Int64, 10)
println("Int64 100x100")
bench_sparse_type(Int64, 100)
println("Int64 1000x1000")
bench_sparse_type(Int64, 1000)

println("Int32 10x10")
bench_sparse_type(Int32, 10)
println("Int32 100x100")
bench_sparse_type(Int32, 100)
println("Int32 1000x1000")
bench_sparse_type(Int32, 1000)

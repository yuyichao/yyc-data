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
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    @inbounds for col in axes(A,2)
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            mul!(C_col, X_col, Aiα, true, true)
        end
    end
    C
end

function spmul4(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β_one = isone(β)
    @inbounds for col in axes(A,2)
        filled = β_one
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            if filled
                mul!(C_col, X_col, Aiα, true, true)
            else
                mul!(C_col, X_col, Aiα, true, β)
            end
            filled = true
        end
        if !filled
            LinearAlgebra._rmul_or_fill!(C_col, β)
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

println("ComplexF64")
println(" 25%")
bench_sparse(zeros(ComplexF64, 100, 100), zeros(ComplexF64, 100, 100),
             sprand(ComplexF64, 100, 100, 0.25))
println(" 5%")
bench_sparse(zeros(ComplexF64, 100, 100), zeros(ComplexF64, 100, 100),
             sprand(ComplexF64, 100, 100, 0.05))
println(" 1%")
bench_sparse(zeros(ComplexF64, 100, 100), zeros(ComplexF64, 100, 100),
             sprand(ComplexF64, 100, 100, 0.01))

println("Float64")
println(" 25%")
bench_sparse(zeros(Float64, 100, 100), zeros(Float64, 100, 100),
             sprand(Float64, 100, 100, 0.25))
println(" 5%")
bench_sparse(zeros(Float64, 100, 100), zeros(Float64, 100, 100),
             sprand(Float64, 100, 100, 0.05))
println(" 1%")
bench_sparse(zeros(Float64, 100, 100), zeros(Float64, 100, 100),
             sprand(Float64, 100, 100, 0.01))

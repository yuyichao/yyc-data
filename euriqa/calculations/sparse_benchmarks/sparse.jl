#!/usr/bin/julia

using LinearAlgebra
using SparseArrays
using SparseArrays: DenseMatrixUnion, SparseMatrixCSCUnion2

using BenchmarkTools

@inline _fast_mul(a, b) = a * b
@inline _fast_mul(a::Union{ComplexF16,ComplexF32,ComplexF64},
                  b::Union{ComplexF16,ComplexF32,ComplexF64}) = muladd(a, b, false)

function spmul_orig(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
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

function spmul_muladd(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    @inbounds for col in axes(A,2), k in nzrange(A, col)
        Aiα = _fast_mul(nzv[k], α)
        rvk = rv[k]
        @simd for multivec_row in axes(X,1)
            C[multivec_row, col] = muladd(X[multivec_row, rvk], Aiα,
                                          C[multivec_row, col])
        end
    end
    C
end

function spmul_view(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    @inbounds for col in axes(A,2)
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = _fast_mul(nzv[k], α)
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            @simd for multivec_row in axes(X,1)
                C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
                                             C_col[multivec_row])
            end
        end
    end
    C
end

function spmul_split(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    rv = rowvals(A)
    nzv = nonzeros(A)
    β_one = isone(β)
    β_zero = iszero(β)
    @inbounds for col in axes(A,2)
        filled = β_one
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = _fast_mul(nzv[k], α)
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            @simd for multivec_row in axes(X,1)
                C_col[multivec_row] = if filled
                    muladd(X_col[multivec_row], Aiα, C_col[multivec_row])
                elseif β_zero
                    _fast_mul(X_col[multivec_row], Aiα)
                else
                    muladd(X_col[multivec_row], Aiα, _fast_mul(C_col[multivec_row], β))
                end
            end
            filled = true
        end
        if !filled
            if β_zero
                fill!(C_col, zero(eltype(C)))
            else
                @simd for multivec_row in axes(X,1)
                    C_col[multivec_row] = _fast_mul(C_col[multivec_row], β)
                end
            end
        end
    end
    C
end

function bench_sparse(C, X, A)
    println("  false")
    @btime spmul_orig($C, $X, $A, true, false)
    @btime spmul_muladd($C, $X, $A, true, false)
    @btime spmul_view($C, $X, $A, true, false)
    @btime spmul_split($C, $X, $A, true, false)

    println("  true")
    @btime spmul_orig($C, $X, $A, true, true)
    @btime spmul_muladd($C, $X, $A, true, true)
    @btime spmul_view($C, $X, $A, true, true)
    @btime spmul_split($C, $X, $A, true, true)

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

println("Float64 10x10")
bench_sparse_type(Float64, 10)
println("Float64 100x100")
bench_sparse_type(Float64, 100)

println("BigFloat 10x10")
bench_sparse_type(BigFloat, 10)
println("BigFloat 100x100")
bench_sparse_type(BigFloat, 100)

println("Int64 10x10")
bench_sparse_type(Int64, 10)
println("Int64 100x100")
bench_sparse_type(Int64, 100)

println("Int32 10x10")
bench_sparse_type(Int32, 10)
println("Int32 100x100")
bench_sparse_type(Int32, 100)

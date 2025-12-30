#!/usr/bin/julia

using LinearAlgebra
using Static
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

function spmul_view(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    rv = rowvals(A)
    nzv = nonzeros(A)
    Xaxes1 = axes(X, 1)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    @inbounds for col in axes(A,2)
        C_col = @view(C[:, col])
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            X_col = @view(X[:, rvk])
            @simd for multivec_row in Xaxes1
                C_col[multivec_row] += X_col[multivec_row] * Aiα
            end
        end
    end
    C
end

@inline _wrapper(A, nrow) = A
struct _Wrapper{T}
    ref::T
    nrow::Int
    offset::Int
end
@inline _wrapper(A::Matrix, nrow) = _Wrapper(A.ref, nrow, 0)
@inline Base.getindex(A::_Wrapper, i) =
    @inbounds Core.memoryrefnew(A.ref, A.offset + i, false)[]
@inline Base.setindex!(A::_Wrapper, v, i) =
    @inbounds Core.memoryrefnew(A.ref, A.offset + i, false)[] = v

@inline _col_view(A, col) = @inbounds @view(A[:, col])
@inline function _col_view(A::_Wrapper, col)
    @inbounds _Wrapper(A.ref, A.nrow, A.offset + (col - 1) * A.nrow)
end

function spmul_view2(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    rv = rowvals(A)
    nzv = nonzeros(A)
    Xaxes1 = axes(X, 1)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    _C = _wrapper(C, mX)
    X = _wrapper(X, mX)
    @inbounds for col in axes(A,2)
        C_col = _col_view(_C, col)
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            X_col = _col_view(X, rvk)
            @simd for multivec_row in Xaxes1
                C_col[multivec_row] += X_col[multivec_row] * Aiα
            end
        end
    end
    C
end

@inline _wrapper2(A, nrow) = A
struct _Wrapper2{T}
    ref::T
    nrow::Int
end
@inline _wrapper2(A::Matrix, nrow) = _Wrapper2(A.ref, nrow)
@inline Base.getindex(A::_Wrapper2, i, j) =
    @inbounds Core.memoryrefnew(A.ref, A.nrow * (j - 1) + i, false)[]
@inline Base.setindex!(A::_Wrapper2, v, i, j) =
    @inbounds Core.memoryrefnew(A.ref, A.nrow * (j - 1) + i, false)[] = v

function spmul_view3(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    rv = rowvals(A)
    nzv = nonzeros(A)
    Xaxes1 = axes(X, 1)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    _C = _wrapper2(C, mX)
    X = _wrapper2(X, mX)
    @inbounds for col in axes(A,2)
        for k in nzrange(A, col)
            Aiα = nzv[k] * α
            rvk = rv[k]
            @simd for multivec_row in Xaxes1
                _C[multivec_row, col] += X[multivec_row, rvk] * Aiα
            end
        end
    end
    C
end

# function spmul_muladd(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
#     mX, nX = size(X)
#     rv = rowvals(A)
#     nzv = nonzeros(A)
#     Xaxes1 = axes(X, 1)
#     β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
#     _C = _wrapper(C, mX)
#     X = _wrapper(X, mX)
#     @inbounds for col in axes(A,2)
#         C_col = _col_view(_C, col)
#         for k in nzrange(A, col)
#             Aiα = nzv[k] * α
#             rvk = rv[k]
#             X_col = _col_view(X, rvk)
#             @simd for multivec_row in Xaxes1
#                 C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
#                                              C_col[multivec_row])
#             end
#         end
#     end
#     C
# end

# function spmul_split(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
#     rv = rowvals(A)
#     nzv = nonzeros(A)
#     Xaxes1 = axes(X, 1)
#     β_one = isone(β)
#     β_zero = iszero(β)
#     if β_zero
#         C_zero = zero(eltype(C))
#     end
#     @inbounds for col in axes(A,2)
#         filled = β_one
#         C_col = @view(C[:, col])
#         for k in nzrange(A, col)
#             Aiα = _fast_mul(nzv[k], α)
#             rvk = rv[k]
#             X_col = @view(X[:, rvk])
#             if filled
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
#                                                  C_col[multivec_row])
#                 end
#             elseif β_zero
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = _fast_mul(X_col[multivec_row], Aiα)
#                 end
#             else
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
#                                                  _fast_mul(C_col[multivec_row], β))
#                 end
#             end
#             filled = true
#         end
#         if !filled
#             if β_zero
#                 fill!(C_col, C_zero)
#             else
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = _fast_mul(C_col[multivec_row], β)
#                 end
#             end
#         end
#     end
#     C
# end

# @inline function _spmul_split(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number,
#                               ::Val{α_one}, ::Val{β_zero}, ::Val{β_one}) where {α_one,β_zero,β_one}
#     rv = rowvals(A)
#     nzv = nonzeros(A)
#     Xaxes1 = axes(X, 1)
#     if β_zero
#         C_zero = zero(eltype(C))
#     end
#     @inbounds for col in axes(A,2)
#         filled = β_one
#         C_col = @view(C[:, col])
#         for k in nzrange(A, col)
#             Aiα = α_one ? nzv[k] : _fast_mul(nzv[k], α)
#             rvk = rv[k]
#             X_col = @view(X[:, rvk])
#             if filled
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
#                                                  C_col[multivec_row])
#                 end
#             elseif β_zero
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = _fast_mul(X_col[multivec_row], Aiα)
#                 end
#             else
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = muladd(X_col[multivec_row], Aiα,
#                                                  _fast_mul(C_col[multivec_row], β))
#                 end
#             end
#             filled = true
#         end
#         if !filled
#             if β_zero
#                 fill!(C_col, C_zero)
#             else
#                 @simd for multivec_row in Xaxes1
#                     C_col[multivec_row] = _fast_mul(C_col[multivec_row], β)
#                 end
#             end
#         end
#     end
# end

# function spmul_split2(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2,
#                       α::Number, β::Number)
#     α_one = isone(α)
#     if isone(β)
#         α_one ? _spmul_split(C, X, A, α, true, Val(true), Val(false), Val(true)) :
#             _spmul_split(C, X, A, α, true, Val(false), Val(false), Val(true))
#     elseif iszero(β)
#         α_one ? _spmul_split(C, X, A, α, false, Val(true), Val(true), Val(false)) :
#             _spmul_split(C, X, A, α, false, Val(false), Val(true), Val(false))
#     else
#         α_one ? _spmul_split(C, X, A, α, β, Val(true), Val(false), Val(false)) :
#             _spmul_split(C, X, A, α, β, Val(false), Val(false), Val(false))
#     end
#     return C
# end

function bench_sparse(C, X, A)
    println("  false")
    @btime spmul_orig($C, $X, $A, true, false)
    @btime spmul_view($C, $X, $A, true, false)
    @btime spmul_view2($C, $X, $A, true, false)
    @btime spmul_view3($C, $X, $A, true, false)
    # @btime spmul_muladd($C, $X, $A, true, false)
    # @btime spmul_split($C, $X, $A, true, false)
    # @btime spmul_split2($C, $X, $A, true, false)

    println("  true")
    @btime spmul_orig($C, $X, $A, true, true)
    @btime spmul_view($C, $X, $A, true, true)
    @btime spmul_view2($C, $X, $A, true, true)
    @btime spmul_view3($C, $X, $A, true, true)
    # @btime spmul_muladd($C, $X, $A, true, true)
    # @btime spmul_split($C, $X, $A, true, true)
    # @btime spmul_split2($C, $X, $A, true, true)

    # println("  static(false)")
    # @btime spmul_orig($C, $X, $A, static(true), static(false))
    # @btime spmul_view($C, $X, $A, static(true), static(false))
    # @btime spmul_view2($C, $X, $A, static(true), static(false))
    # @btime spmul_view3($C, $X, $A, static(true), static(false))
    # # @btime spmul_muladd($C, $X, $A, static(true), static(false))
    # # @btime spmul_split($C, $X, $A, static(true), static(false))
    # # @btime spmul_split2($C, $X, $A, static(true), static(false))

    # println("  static(true)")
    # @btime spmul_orig($C, $X, $A, static(true), static(true))
    # @btime spmul_view($C, $X, $A, static(true), static(true))
    # @btime spmul_view2($C, $X, $A, static(true), static(true))
    # @btime spmul_view3($C, $X, $A, static(true), static(true))
    # # @btime spmul_muladd($C, $X, $A, static(true), static(true))
    # # @btime spmul_split($C, $X, $A, static(true), static(true))
    # # @btime spmul_split2($C, $X, $A, static(true), static(true))

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
    println(" 0%")
    bench_sparse(zeros(ElType, sz, sz), zeros(ElType, sz, sz),
                 spzeros(ElType, sz, sz))
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

#!/usr/bin/julia

using LinearAlgebra
using SparseArrays
using SparseArrays: DenseMatrixUnion, SparseMatrixCSCUnion2

@inline _fast_mul(a, b) = a * b
@inline _fast_mul(a::Union{ComplexF16,ComplexF32,ComplexF64},
                  b::Union{ComplexF16,ComplexF32,ComplexF64}) = muladd(a, b, false)

function spmul_orig(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    mA, nA = size(A)
    mC, nC = size(C)
    nX == mA ||
        throw(DimensionMismatch("second dimension of X, $nX, does not match the first dimension of A, $(mA)"))
    mX == mC ||
        throw(DimensionMismatch("first dimension of X, $mX, does not match the first dimension of C, $(mC)"))
    nA == nC ||
        throw(DimensionMismatch("second dimension of A, $(nA), does not match the second dimension of C, $(nC)"))
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

@inline _wrap_matrix(A, nrow) = A
struct _ImmutableMatrix{T}
    ref::T
    nrow::Int
end
@inline _wrap_matrix(A::Matrix, nrow) = _ImmutableMatrix(A.ref, nrow)
@inline Base.getindex(A::_ImmutableMatrix, i, j) =
    @inbounds Core.memoryrefnew(A.ref, A.nrow * (j - 1) + i, false)[]
@inline Base.setindex!(A::_ImmutableMatrix, v, i, j) =
    @inbounds Core.memoryrefnew(A.ref, A.nrow * (j - 1) + i, false)[] = v

function spmul_view(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    Xaxes1 = axes(X, 1)
    mA, nA = size(A)
    Aaxes2 = axes(A, 2)
    mC, nC = size(C)
    nX == mA ||
        throw(DimensionMismatch("second dimension of X, $nX, does not match the first dimension of A, $(mA)"))
    mX == mC ||
        throw(DimensionMismatch("first dimension of X, $mX, does not match the first dimension of C, $(mC)"))
    nA == nC ||
        throw(DimensionMismatch("second dimension of A, $(nA), does not match the second dimension of C, $(nC)"))
    rv = rowvals(A)
    nzv = nonzeros(A)
    β != one(β) && LinearAlgebra._rmul_or_fill!(C, β)
    _C = _wrap_matrix(C, mX)
    X = _wrap_matrix(X, mX)
    @inbounds for col in Aaxes2
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

function spmul_muladd(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    mX, nX = size(X)
    Xaxes1 = axes(X, 1)
    mA, nA = size(A)
    Aaxes2 = axes(A, 2)
    mC, nC = size(C)
    nX == mA ||
        throw(DimensionMismatch("second dimension of X, $nX, does not match the first dimension of A, $(mA)"))
    mX == mC ||
        throw(DimensionMismatch("first dimension of X, $mX, does not match the first dimension of C, $(mC)"))
    nA == nC ||
        throw(DimensionMismatch("second dimension of A, $(nA), does not match the second dimension of C, $(nC)"))
    rv = rowvals(A)
    nzv = nonzeros(A)
    isone(β) || LinearAlgebra._rmul_or_fill!(C, β)
    if α isa Bool && !α
        return C
    end
    _C = _wrap_matrix(C, mX)
    X = _wrap_matrix(X, mX)
    @inbounds for col in Aaxes2
        for k in nzrange(A, col)
            Aiα = α isa Bool ? nzv[k] : nzv[k] * α
            rvk = rv[k]
            @simd for multivec_row in Xaxes1
                _C[multivec_row, col] = muladd(X[multivec_row, rvk], Aiα,
                                               _C[multivec_row, col])
            end
        end
    end
    C
end

@noinline function _spmul_error(mX, nX, mA, nA, mC, nC)
    nX == mA ||
        throw(DimensionMismatch("second dimension of X, $nX, does not match the first dimension of A, $(mA)"))
    mX == mC ||
        throw(DimensionMismatch("first dimension of X, $mX, does not match the first dimension of C, $(mC)"))
    nA == nC ||
        throw(DimensionMismatch("second dimension of A, $(nA), does not match the second dimension of C, $(nC)"))
end

function _spmul_split(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number, Annz, ::Val{β_zero}, ::Val{β_one}, ::Val{Prefill}) where {β_zero,β_one,Prefill}
    mX, nX = size(X)
    Xaxes1 = axes(X, 1)
    mA, nA = size(A)
    Aaxes2 = axes(A, 2)
    mC, nC = size(C)
    if !(nX == mA & mX == mC & nA == nC)
        _spmul_error(mX, nX, mA, nA, mC, nC)
    end
    rv = rowvals(A)
    nzv = nonzeros(A)
    if β isa Bool
        β = β_one
    end
    empty_mul = (α isa Bool && !α) | Annz == 0

    if Prefill || empty_mul
        β_one || LinearAlgebra._rmul_or_fill!(C, β)
        if empty_mul
            return C
        end
    end
    if !Prefill && β_zero
        C_zero = zero(eltype(C))
    end

    _C = _wrap_matrix(C, mX)
    X = _wrap_matrix(X, mX)
    @inbounds for col in Aaxes2
        nzrng = nzrange(A, col)
        if isempty(nzrng)
            if Prefill
                # Already filled
            elseif β_zero
                @simd for multivec_row in Xaxes1
                    _C[multivec_row, col] = C_zero
                end
            else
                @simd for multivec_row in Xaxes1
                    _C[multivec_row, col] = _fast_mul(_C[multivec_row, col], β)
                end
            end
            continue
        end
        filled = Prefill
        for k in nzrng
            Aiα = α isa Bool ? nzv[k] : nzv[k] * α
            rvk = rv[k]
            if !filled
                if β_zero
                    @simd for multivec_row in Xaxes1
                        _C[multivec_row, col] = _fast_mul(X[multivec_row, rvk], Aiα)
                    end
                else
                    @simd for multivec_row in Xaxes1
                        _C[multivec_row, col] = muladd(X[multivec_row, rvk], Aiα,
                                                       _C[multivec_row, col] * β)
                    end
                end
                filled = true
            else
                @simd for multivec_row in Xaxes1
                    _C[multivec_row, col] = muladd(X[multivec_row, rvk], Aiα,
                                                   _C[multivec_row, col])
                end
            end
        end
    end
    C
end

function spmul_split(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2,
                     α::Number, β::Number)
    Annz = nnz(A)
    if isone(β)
        _spmul_split(C, X, A, α, true, Annz, Val(false), Val(true), Val(true))
    elseif iszero(β)
        (isbitstype(eltype(C)) && Annz <= 10) ?
            _spmul_split(C, X, A, α, false, Annz, Val(true), Val(false), Val(true)) :
            _spmul_split(C, X, A, α, false, Annz, Val(true), Val(false), Val(false))
    else
        (isbitstype(eltype(C)) && Annz <= 10) ?
            _spmul_split(C, X, A, α, β, Annz, Val(false), Val(false), Val(true)) :
            _spmul_split(C, X, A, α, β, Annz, Val(false), Val(false), Val(false))
    end
    return C
end

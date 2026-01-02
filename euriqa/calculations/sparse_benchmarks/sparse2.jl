#!/usr/bin/julia

using LinearAlgebra
using SparseArrays
using SparseArrays: AdjOrTrans, DenseMatrixUnion, SparseMatrixCSCUnion2

_fix_size(M, nrow, ncol) = M
_fix_size(M) = _fix_size(M, size(M)...)

@static if VERSION >= v"1.11"
    # An immutable fixed size wrapper for matrices to work around
    # the performance issue caused by https://github.com/JuliaLang/julia/issues/60409
    struct _FixedSizeMatrix{Trans,R}
        ref::R
        nrow::Int
        ncol::Int
        function _FixedSizeMatrix{Trans}(ref::R, nrow, ncol) where {Trans,R}
            new{Trans,R}(ref, nrow, ncol)
        end
    end
    @inline Base.getindex(A::_FixedSizeMatrix{'N'}, i, j) =
        @inbounds Core.memoryrefnew(A.ref, A.nrow * (j - 1) + i, false)[]
    @inline Base.setindex!(A::_FixedSizeMatrix{'N'}, v, i, j) =
        @inbounds Core.memoryrefnew(A.ref, A.nrow * (j - 1) + i, false)[] = v

    @inline Base.getindex(A::_FixedSizeMatrix{'T'}, i, j) =
        @inbounds transpose(Core.memoryrefnew(A.ref, A.ncol * (i - 1) + j, false)[])
    @inline Base.setindex!(A::_FixedSizeMatrix{'T'}, v, i, j) =
        @inbounds Core.memoryrefnew(A.ref, A.ncol * (i - 1) + j, false)[] = transpose(v)

    @inline Base.getindex(A::_FixedSizeMatrix{'C'}, i, j) =
        @inbounds adjoint(Core.memoryrefnew(A.ref, A.ncol * (i - 1) + j, false)[])
    @inline Base.setindex!(A::_FixedSizeMatrix{'C'}, v, i, j) =
        @inbounds Core.memoryrefnew(A.ref, A.ncol * (i - 1) + j, false)[] = adjoint(v)

    @inline _fix_size(A::Matrix, nrow, ncol) = _FixedSizeMatrix{'N'}(A.ref, nrow, ncol)
    @inline _fix_size(A::Transpose{<:Any,<:Matrix}, nrow, ncol) =
        _FixedSizeMatrix{'T'}(A.parent.ref, nrow, ncol)
    @inline _fix_size(A::Adjoint{<:Any,<:Matrix}, nrow, ncol) =
        _FixedSizeMatrix{'C'}(A.parent.ref, nrow, ncol)
end

# Slow non-inlined functions for throwing the error without messing up the caller
@noinline function _matmul_size_error(mC, nC, mA, nA, mB, nB, At, Bt)
    if At == 'N'
        Anames = "first", "second"
    else
        Anames = "second", "first"
    end
    if Bt == 'N'
        Bnames = "first", "second"
    else
        Bnames = "second", "first"
    end
    nA == mB ||
        throw(DimensionMismatch("$(Anames[2]) dimension of A, $nA, does not match the $(Bnames[1]) dimension of B, $mB"))
    mA == mC ||
        throw(DimensionMismatch("$(Anames[1]) dimension of A, $mA, does not match the first dimension of C, $mC"))
    nB == nC ||
        throw(DimensionMismatch("$(Bnames[2]) dimension of B, $nB, does not match the second dimension of C, $nC"))
    # unreachable
    throw(DimensionMismatch("Unknown dimension mismatch"))
end

@inline function _matmul_size(C, A, B, ::Val{At}, ::Val{Bt}) where {At,Bt}
    mC = size(C, 1)
    nC = size(C, 2)
    mA = size(A, 1)
    nA = size(A, 2)
    mB = size(B, 1)
    nB = size(B, 2)

    _mA, _nA = At == 'N' ? (mA, nA) : (nA, mA)
    _mB, _nB = Bt == 'N' ? (mB, nB) : (nB, mB)

    if (_nA != _mB) | (_mA != mC) | (_nB != nC)
        _matmul_size_error(mC, nC, _mA, _nA, _mB, _nB, At, Bt)
    end
    return mC, nC, mA, nA, mB, nB
end

@inline _matmul_size_AB(C, A, B) = _matmul_size(C, A, B, Val('N'), Val('N'))
@inline _matmul_size_AtB(C, A, B) = _matmul_size(C, A, B, Val('T'), Val('N'))
@inline _matmul_size_ABt(C, A, B) = _matmul_size(C, A, B, Val('N'), Val('T'))

function spmul_orig(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    Aax2 = axes(A, 2)
    Xax1 = axes(X, 1)
    mC, nC, mX, nX, mA, nA = _matmul_size_AB(C, X, A)
    rv = rowvals(A)
    nzv = nonzeros(A)
    isone(β) || LinearAlgebra._rmul_or_fill!(C, β)
    if α isa Bool && !α
        return
    end
    C = _fix_size(C, mC, nC)
    X = _fix_size(X, mX, nX)
    @inbounds for col in Aax2, k in nzrange(A, col)
        Aiα = α isa Bool ? nzv[k] : nzv[k] * α
        rvk = rv[k]
        @simd for multivec_row in Xax1
            C[multivec_row, col] = muladd(X[multivec_row, rvk], Aiα,
                                          C[multivec_row, col])
        end
    end
end

function spmul_adj_orig(C::StridedMatrix, X::AdjOrTrans{<:Any,<:DenseMatrixUnion}, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    Xax1 = axes(X, 1)
    Cax2 = axes(C, 2)
    mC, nC, mX, nX, mA, nA = _matmul_size_AB(C, X, A)
    rv = rowvals(A)
    nzv = nonzeros(A)
    isone(β) || LinearAlgebra._rmul_or_fill!(C, β)
    if α isa Bool && !α
        return
    end
    C = _fix_size(C, mC, nC)
    X = _fix_size(X, mX, nX)
    @inbounds for multivec_row in Xax1, col in Cax2
        nzrng = nzrange(A, col)
        if isempty(nzrng)
            continue
        end
        tmp = C[multivec_row, col]
        for k in nzrng
            tmp = muladd(X[multivec_row, rv[k]],
                         (α isa Bool ? nzv[k] : nzv[k] * α), tmp)
        end
        C[multivec_row, col] = tmp
    end
end

@inline _fast_mul(a, b) = a * b
@inline _fast_mul(a::Union{ComplexF16,ComplexF32,ComplexF64},
                  b::Union{ComplexF16,ComplexF32,ComplexF64}) = muladd(a, b, false)

@inline _col_fill!(C, ccol, α, ax) = @inbounds @simd for i in ax
    C[i, ccol] = α
end
@inline _col_mul!(C, ccol, X, xcol, α, ax) = @inbounds @simd for i in ax
    C[i, ccol] = _fast_mul(X[i, xcol], α)
end
@inline _col_muladd!(C, ccol, X, xcol, α, ax) = @inbounds @simd for i in ax
    C[i, ccol] = muladd(X[i, xcol], α, C[i, ccol])
end
@inline _col_muladdmul!(C, ccol, X, xcol, α, β, ax) = @inbounds @simd for i in ax
    C[i, ccol] = muladd(X[i, xcol], α, _fast_mul(C[i, ccol], β))
end

function spmul_split(C::StridedMatrix, X::DenseMatrixUnion, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    Aax2 = axes(A, 2)
    Xax1 = axes(X, 1)
    mC, nC, mX, nX, mA, nA = _matmul_size_AB(C, X, A)
    rv = rowvals(A)
    nzv = nonzeros(A)

    _C = _fix_size(C, mC, nC)
    X = _fix_size(X, mX, nX)

    αBool = α isa Bool
    ElType = eltype(C)

    if isone(β)
        @goto mul_only
    end
    if (isbitstype(ElType) && length(rv) < 10) || α === false
        LinearAlgebra._rmul_or_fill!(C, β)
        @label mul_only
        if α !== false
            @inbounds for col in Aax2, k in nzrange(A, col)
                _col_muladd!(_C, col, X, rv[k], αBool ? nzv[k] : nzv[k] * α, Xax1)
            end
        end
        return
    end

    if iszero(β)
        C_zero = zero(ElType)
        @inbounds for col in Aax2
            nzrng = nzrange(A, col)
            if isempty(nzrng)
                _col_fill!(_C, col, C_zero, Xax1)
                continue
            end
            filled = false
            for k in nzrng
                Aiα = αBool ? nzv[k] : nzv[k] * α
                rvk = rv[k]
                if !filled
                    _col_mul!(_C, col, X, rvk, Aiα, Xax1)
                    filled = true
                else
                    _col_muladd!(_C, col, X, rvk, Aiα, Xax1)
                end
            end
        end
    else
        @inbounds for col in Aax2
            nzrng = nzrange(A, col)
            if isempty(nzrng)
                _col_mul!(_C, col, _C, col, β, Xax1)
                continue
            end
            filled = false
            for k in nzrng
                Aiα = αBool ? nzv[k] : nzv[k] * α
                rvk = rv[k]
                if !filled
                    _col_muladdmul!(_C, col, X, rvk, Aiα, β, Xax1)
                    filled = true
                else
                    _col_muladd!(_C, col, X, rvk, Aiα, Xax1)
                end
            end
        end
    end
end

@inline function _spmul_adj_accum(X, i, rv, nzv, nzrng, α, accum)
    @inbounds @simd for k in nzrng
        accum = muladd(X[i, rv[k]], (α isa Bool ? nzv[k] : nzv[k] * α), accum)
    end
    return accum
end

function spmul_adj_order(C::StridedMatrix, X::AdjOrTrans{<:Any,<:DenseMatrixUnion}, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    Xax1 = axes(X, 1)
    Cax2 = axes(C, 2)
    mC, nC, mX, nX, mA, nA = _matmul_size_AB(C, X, A)
    rv = rowvals(A)
    nzv = nonzeros(A)
    isone(β) || LinearAlgebra._rmul_or_fill!(C, β)
    C = _fix_size(C, mC, nC)
    X = _fix_size(X, mX, nX)
    if (α isa Bool && !α) || isempty(Xax1)
        return
    end
    @inbounds for col in Cax2
        nzrng = nzrange(A, col)
        if isempty(nzrng)
            continue
        end
        for multivec_row in Xax1
            tmp = C[multivec_row, col]
            for k in nzrng
                tmp = muladd(X[multivec_row, rv[k]],
                             (α isa Bool ? nzv[k] : nzv[k] * α), tmp)
            end
            C[multivec_row, col] = tmp
        end
    end
end

function spmul_adj_split(C::StridedMatrix, X::AdjOrTrans{<:Any,<:DenseMatrixUnion}, A::SparseMatrixCSCUnion2, α::Number, β::Number)
    Xax1 = axes(X, 1)
    Cax2 = axes(C, 2)
    mC, nC, mX, nX, mA, nA = _matmul_size_AB(C, X, A)
    rv = rowvals(A)
    nzv = nonzeros(A)

    ElType = eltype(C)
    _C = _fix_size(C, mC, nC)
    X = _fix_size(X, mX, nX)

    @inbounds if isone(β)
        if α isa Bool && !α
            return
        end
        for col in Cax2
            nzrng = nzrange(A, col)
            if isempty(nzrng)
                continue
            end
            for multivec_row in Xax1
                _C[multivec_row, col] = _spmul_adj_accum(X, multivec_row, rv, nzv, nzrng,
                                                         α, _C[multivec_row, col])
            end
        end
    else
        if α isa Bool && !α
            LinearAlgebra._rmul_or_fill!(C, β)
            return
        end
        if iszero(β)
            C_zero = zero(ElType)
            for col in Cax2
                nzrng = nzrange(A, col)
                for multivec_row in Xax1
                    _C[multivec_row, col] = _spmul_adj_accum(X, multivec_row, rv, nzv,
                                                             nzrng, α, C_zero)
                end
            end
        else
            for col in Cax2
                nzrng = nzrange(A, col)
                for multivec_row in Xax1
                    _C[multivec_row, col] = _spmul_adj_accum(X, multivec_row, rv, nzv,
                                                             nzrng, α,
                                                             _C[multivec_row, col] * β)
                end
            end
        end
    end
end

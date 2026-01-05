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

@noinline function spmul_adj_orig2(C, X, A, colptr, sz)
    nzv = A.nzval
    @inbounds begin
        for i in Base.unchecked_oneto(sz)
            nzrng = Base.unchecked_oneto(colptr[])
            tmp = C[i]
            for k in nzrng
                tmp = muladd(X[k], nzv[k], tmp)
            end
            C[i] = tmp
        end
    end
end

@noinline function spmul_adj_order(C, X, A, colptr, sz)
    nzv = A.nzval
    @inbounds begin
        nzrng = Base.unchecked_oneto(colptr[])
        for i in Base.unchecked_oneto(sz)
            tmp = C[i]
            for k in nzrng
                tmp = muladd(X[k], nzv[k], tmp)
            end
            C[i] = tmp
        end
    end
end

using BenchmarkTools
using StaticArrays
using InteractiveUtils

function save_code(f, args, prefix)
    types = Base.typesof(args...)
    open(io->code_llvm(io, f, types, raw=true, dump_module=true,
                       optimize=false, debuginfo=:none),
         "$(prefix)-raw_noopt.ll", "w")
    open(io->code_llvm(io, f, types, dump_module=true, optimize=false, debuginfo=:none),
         "$(prefix)-noopt.ll", "w")
    open(io->code_llvm(io, f, types, dump_module=true, optimize=true, debuginfo=:none),
         "$(prefix)-opt.ll", "w")
    open(io->code_native(io, f, types, dump_module=true, debuginfo=:none),
         "$(prefix)-opt.s", "w")
end

function bench_sparse_type(ElType, sz)
    C = (zeros(ElType, sz))
    X = MVector{sz}(zeros(ElType, sz))
    A = (;sz = sz, nzval = (zeros(ElType, sz)))
    @btime spmul_adj_orig2($C, $X, $A, $(Ref(sz)), $sz)
    @btime spmul_adj_order($C, $X, $A, $(Ref(sz)), $sz)
    save_code(spmul_adj_orig2, (C, X, A, Ref(sz), sz), "orig2")
    save_code(spmul_adj_order, (C, X, A, Ref(sz), sz), "order")
    return
end

bench_sparse_type(ComplexF64, 300)

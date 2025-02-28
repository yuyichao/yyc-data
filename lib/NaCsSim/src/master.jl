#!/usr/bin/julia

module Master

using LinearAlgebra
using SparseArrays

import OrdinaryDiffEq, DiffEqCallbacks

function foreach_nz(cb, m::SparseMatrixCSC)
    for col in 1:size(m, 2)
        for r in nzrange(m, col)
            row = rowvals(m)[r]
            v = nonzeros(m)[r]
            cb(row, col, v)
        end
    end
end

@generated function map_op!(::Val{Map}, tgt, src) where Map
    ex = quote
    end

    max_iin = 0
    max_iout = 0
    for (iout, iin, v) in Map
        max_iin = max(max_iin, iin)
        max_iout = max(max_iout, iout)
    end

    num_map = length(Map)
    minpos_in = fill(num_map, max_iin)
    minpos_out = fill(num_map, max_iout)
    maxpos_out = fill(1, max_iout)

    for i in 1:num_map
        iout, iin, v = Map[i]
        minpos_in[iin] = min(minpos_in[iin], i)
        minpos_out[iout] = min(minpos_out[iout], i)
        maxpos_out[iout] = max(maxpos_out[iout], i)
    end

    invars = [gensym() for i in 1:max_iin]
    outvars = [gensym() for i in 1:max_iout]

    for i in 1:num_map
        iout, iin, v = Map[i]
        invar = invars[iin]
        outvar = outvars[iout]
        if i == minpos_in[iin]
            push!(ex.args, :($invar = @inbounds src[$iin]))
        end
        if i == minpos_out[iout]
            push!(ex.args, :($outvar = @inbounds tgt[$iout]))
        end
        push!(ex.args, :($outvar = muladd($invar, $v, $outvar)))
        if i == maxpos_out[iout]
            push!(ex.args, :(@inbounds tgt[$iout] = $outvar))
        end
    end
    push!(ex.args, :(return))
    return ex
end

struct System{T<:Real,N,D,JMap}
    op_dagger::SparseMatrixCSC{Complex{T},Int}
    data::D
    function System{T,_N}(Js, data::D, N=nothing) where {T<:Real, _N, D}
        if _N !== nothing
            N = _N
        end
        linidx = LinearIndices((N, N))
        jmap = Dict{NTuple{2,Int},T}()
        function add_jterm(xi, yi, xo, yo, v)
            key = linidx[xi, yi], linidx[xo, yo]
            if haskey(jmap, key)
                jmap[key] += v
            else
                jmap[key] = v
            end
            return
        end
        for j in Js
            foreach_nz(j) do x1, y1, v1
                foreach_nz(j) do x2, y2, v2
                    add_jterm(y1, y2, x1, x2, v1 * v2)
                end
            end
        end
        jop = [(iout, iin, v) for ((iin, iout), v) in jmap]
        # Group the ones with the same output together
        sort!(jop)
        return new{T,_N,D,tuple(jop...)}(spzeros(T, N, N), data)
    end
end

struct SystemCoherent{OP<:AbstractMatrix,N,D}
    op::OP
    data::D
    function SystemCoherent{OP,_N}(data::D, N=nothing) where {OP,_N,D}
        if _N !== nothing
            N = _N
        end
        return new{OP,_N,D}(zeros(N, N), data)
    end
end

function init! end
function update! end

@inline function get_nrow(M, ::Val{nothing})
    return size(M, 1)
end
@inline function get_nrow(M, ::Val{nrow}) where nrow
    return nrow
end

@inline function do_mul!(result::DenseMatrix, B::DenseMatrix,
                         M::SparseMatrixCSC, ::Val{_nrow}) where _nrow
    prev_colptr = 1
    nrow = get_nrow(M, Val(_nrow))
    @inbounds for col in 1:nrow
        filled = false
        colptr = M.colptr[col + 1]
        # Manually compute the linear index.
        # Somehow the invariance of the matrix size doesn't work anymore
        # and the dimension of the matrix is being loaded in every loop iterations
        col_offset = (col - 1) * nrow
        for i in prev_colptr:colptr - 1
            val = M.nzval[i]
            row = M.rowval[i]
            row_offset = (row - 1) * nrow
            if filled
                @simd ivdep for j in 1:nrow
                    result[j + col_offset] =
                        muladd(val, B[j + row_offset], result[j + col_offset])
                end
            else
                @simd ivdep for j in 1:nrow
                    result[j + col_offset] = val * B[j + row_offset]
                end
            end
            filled = true
        end
        prev_colptr = colptr
        if !filled
            @simd ivdep for j in 1:nrow
                result[j + col_offset] = 0
            end
        end
    end
end

@inline function dmaster(t, rho_data::DenseMatrix, drho_data::DenseMatrix,
                         sys::System{T,_nrow,D,JMap}, drive) where {T,_nrow,D,JMap}
    @inline update!(drive, sys, t)
    @inline do_mul!(drho_data, rho_data, sys.op_dagger, Val{_nrow}())
    # compute -i * drho^dagger + i * drho
    nrow = get_nrow(sys.op_dagger, Val(_nrow))
    @inbounds for i in 1:nrow
        drho_data[i + (i - 1) * nrow] = -2 * imag(drho_data[i + (i - 1) * nrow])
        for j in (i + 1):nrow
            idx1 = j + (i - 1) * nrow
            idx2 = i + (j - 1) * nrow
            v1 = drho_data[idx1]
            v2 = drho_data[idx2]

            re_o = -imag(v1) - imag(v2)
            im_o = real(v2) - real(v1)

            drho_data[idx1] = complex(re_o, -im_o)
            drho_data[idx2] = complex(re_o, im_o)
        end
    end
    @inline map_op!(Val{JMap}(), drho_data, rho_data)
    return
end

function _mul_hermitian!(result::DenseVector, M::SparseMatrixCSC,
                         B::DenseVector, ::Val{_nrow}) where _nrow
    nrow = get_nrow(M, Val(_nrow))
    colptrs = M.colptr
    rowval = M.rowval
    nzval = M.nzval
    prev_colptr = 1
    @inbounds for col in 1:nrow
        r = zero(eltype(result))
        colptr = colptrs[col + 1]
        for i in prev_colptr:colptr - 1
            Bv = B[rowval[i]]
            r = muladd(conj(nzval[i]), complex(-imag(Bv), real(Bv)), r)
        end
        result[col] = r
        prev_colptr = colptr
    end
end
_mul_hermitian!(result::DenseVector, M, B::DenseVector, _) = @inline mul!(result, M, B)

@inline function dschroedinger(t, ψ::DenseVector, dψ::DenseVector,
                               sys::SystemCoherent{OP,_nrow}, drive) where {OP,_nrow}
    @inline update!(drive, sys, t)
    @inline _mul_hermitian!(dψ, sys.op, ψ, Val(_nrow))
end

function integrate(tspan, df_::F, x0, fout::FO;
                   alg = OrdinaryDiffEq.DP5(),
                   save_everystep = false, saveat=tspan,
                   callback = nothing, kwargs...) where {F,FO}

    function fout_(x, t, integrator)
        @inline fout(t, x)
    end

    tType = float(eltype(tspan))
    out_type = pure_inference(fout, Tuple{tType,typeof(x0)})

    out = DiffEqCallbacks.SavedValues(tType,out_type)

    scb = DiffEqCallbacks.SavingCallback(fout_, out, saveat=saveat,
                                         save_everystep=save_everystep,
                                         save_start = false,
                                         tdir = one(eltype(tspan)))

    prob = OrdinaryDiffEq.ODEProblem{true}(df_, x0, (convert(tType, tspan[1]),
                                                     convert(tType, tspan[end])))

    full_cb = OrdinaryDiffEq.CallbackSet(callback, scb)

    sol = OrdinaryDiffEq.solve(
                prob,
                alg;
                reltol = 1.0e-6,
                abstol = 1.0e-8,
                save_everystep = false, save_start = false,
                save_end = false,
                callback=full_cb, kwargs...)
    return out.t, out.saveval
end

Base.@pure pure_inference(fout, T) = Core.Compiler.return_type(fout, T)

function evolve(drive, sys::System, x0::DenseMatrix, tspan;
                fout::FO=nothing, kws...) where FO
    init!(drive, sys)
    function dmaster_(dx, x, p, t)
        dmaster(t, x, dx, sys, drive)
    end
    if fout === nothing
        fout = (t, x) -> copy(x)
    end
    return integrate(tspan, dmaster_, x0, fout; kws...)
end

function evolve(drive, sys::SystemCoherent, x0::DenseVector, tspan;
                fout::FO=nothing, kws...) where FO
    init!(drive, sys)
    function dschroedinger_(dx, x, p, t)
        dschroedinger(t, x, dx, sys, drive)
    end
    if fout === nothing
        fout = (t, x) -> copy(x)
    end
    return integrate(tspan, dschroedinger_, x0, fout; kws...)
end

end

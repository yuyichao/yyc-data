#!/usr/bin/julia

include("sparse2.jl")

using BenchmarkTools

function bench_sparse(C, X, A)
    println("  adjoint false")
    @btime spmul_adj_orig($C, $(adjoint(X)), $A, true, false)
    @btime spmul_adj_order($C, $(adjoint(X)), $A, true, false)
    @btime spmul_adj_split($C, $(adjoint(X)), $A, true, false)
    println("  transpose false")
    @btime spmul_adj_orig($C, $(transpose(X)), $A, true, false)
    @btime spmul_adj_order($C, $(transpose(X)), $A, true, false)
    @btime spmul_adj_split($C, $(transpose(X)), $A, true, false)

    println("  adjoint true")
    @btime spmul_adj_orig($C, $(adjoint(X)), $A, true, true)
    @btime spmul_adj_order($C, $(adjoint(X)), $A, true, true)
    @btime spmul_adj_split($C, $(adjoint(X)), $A, true, true)
    println("  transpose true")
    @btime spmul_adj_orig($C, $(transpose(X)), $A, true, true)
    @btime spmul_adj_order($C, $(transpose(X)), $A, true, true)
    @btime spmul_adj_split($C, $(transpose(X)), $A, true, true)
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
    A = sprand(ElType, sz, sz, 0.01)
    while nnz(A) == 0
        A = sprand(ElType, sz, sz, 0.01)
    end
    bench_sparse(zeros(ElType, sz, sz), zeros(ElType, sz, sz), A)
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

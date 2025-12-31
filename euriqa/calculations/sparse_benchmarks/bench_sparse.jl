#!/usr/bin/julia

include("sparse.jl")

using BenchmarkTools
using YAML

function bench_sparse_param(C, X, A, α, β)
    return Dict{String,Any}("orig"=>@btimed(spmul_orig($C, $X, $A, $α, $β)),
                            "view"=>@btimed(spmul_view($C, $X, $A, $α, $β)),
                            "muladd"=>@btimed(spmul_muladd($C, $X, $A, $α, $β)),
                            "split"=>@btimed(spmul_split($C, $X, $A, $α, $β)))
end

function bench_sparse_matrix(C, X, A)
    ET = eltype(C)
    vs = Any[false, true, 0, 1, 2]
    res = Dict{String,Any}[]
    for α in vs
        for β in vs
            b = bench_sparse_param(f!, C, X, A, ET(α), ET(β))
            b["α"] = α
            b["β"] = β
            push!(res, b)
        end
    end
    return res
end

function bench_sparse_type_sz(ElType, sz)
    C = zeros(ElType, sz, sz)
    X = zeros(ElType, sz, sz)

    println(" 25%")
    b25 = bench_sparse_matrix(C, X, sprand(ElType, sz, sz, 0.25))
    println(" 5%")
    b5 = bench_sparse_matrix(C, X, sprand(ElType, sz, sz, 0.05))
    println(" 1%")
    A = sprand(ElType, sz, sz, 0.01)
    while nnz(A) == 0
        A = sprand(ElType, sz, sz, 0.01)
    end
    b1 = bench_sparse_matrix(C, X, A)
    println(" 0%")
    b0 = bench_sparse_matrix(C, X, spzeros(ElType, sz, sz))
    println(" diag")
    bdiag = bench_sparse_matrix(C, X, sparse(1:sz, 1:sz, ones(ElType, sz)))

    return Dict("25"=>b25, "5"=>b5, "1"=>b1, "0"=>b0, "diag"=>bdiag)
end

function bench_sparse_type(ElType)
    return Dict("10"=>bench_sparse_type_sz(ElType, 10),
                "100"=>bench_sparse_type_sz(ElType, 100))
end

function bench_sparse()
    return Dict("ComplexF64"=>bench_sparse_type(ComplexF64),
                "Float64"=>bench_sparse_type(Float64),
                "Int64"=>bench_sparse_type(Int64),
                "BigFloat"=>bench_sparse_type(BigFloat))
end

YAML.write_file(ARGS[1], bench_sparse())

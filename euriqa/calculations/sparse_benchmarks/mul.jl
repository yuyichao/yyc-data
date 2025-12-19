#!/usr/bin/julia

using LinearAlgebra

using BenchmarkTools

@inline mul1(C, A, α) = mul!(C, A, α)
function mul2(C, A, α)
    @inbounds @simd for i in 1:length(C)
        C[i] += A[i] * α
    end
end

function bench_mul(C, A, α)
    @btime mul1($C, $A, $α)
    @btime mul2($C, $A, $α)
    return
end

println("ComplexF64")
println(" 100")
bench_mul(zeros(ComplexF64, 100), zeros(ComplexF64, 100), 0.0im)
println(" 1000")
bench_mul(zeros(ComplexF64, 1000), zeros(ComplexF64, 1000), 0.0im)

println("Float64")
println(" 100")
bench_mul(zeros(Float64, 100), zeros(Float64, 100), 0.0)
println(" 1000")
bench_mul(zeros(Float64, 1000), zeros(Float64, 1000), 0.0)

println("ComplexF64 sub")
println(" 100")
bench_mul(@view(zeros(ComplexF64, 100, 1)[:, 1]),
          @view(zeros(ComplexF64, 100, 1)[:, 1]), 0.0im)
println(" 1000")
bench_mul(@view(zeros(ComplexF64, 1000, 1)[:, 1]),
          @view(zeros(ComplexF64, 1000, 1)[:, 1]), 0.0im)

println("Float64 sub")
println(" 100")
bench_mul(@view(zeros(Float64, 100, 1)[:, 1]),
          @view(zeros(Float64, 100, 1)[:, 1]), 0.0)
println(" 1000")
bench_mul(@view(zeros(Float64, 1000, 1)[:, 1]),
          @view(zeros(Float64, 1000, 1)[:, 1]), 0.0)

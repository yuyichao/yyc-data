#!/usr/bin/julia

include("ind-modulation.jl")

using Test

@testset "Test area" begin
    for _ in 1:1000
        dτ = rand() + 0.2
        ω = 2 + rand()
        ωm = 2.5 + rand()
        Ω1s = rand(20) .* 2 .- 1
        Ω2s = rand(20) .* 2 .- 1

        a11 = compute_area2_am(dτ, dτ, ω, ωm, Ω1s, Ω1s)
        a12 = compute_area2_am(dτ, dτ, ω, ωm, Ω1s, Ω2s)
        a21 = compute_area2_am(dτ, dτ, ω, ωm, Ω2s, Ω1s)
        a22 = compute_area2_am(dτ, dτ, ω, ωm, Ω2s, Ω2s)
        a12_12 = compute_area2_am(dτ, dτ, ω, ωm, Ω1s .+ Ω2s, Ω1s .+ Ω2s)
        @test a11 + a12 + a21 + a22 ≈ a12_12
    end
end

#!/usr/bin/julia

include("util.jl")

using BenchmarkTools

const param = Params{Float64}(1.0, 0.01, 1.0, 110, (24, 24, 24))
const state = State(param)

@btime sweep!(state)
@btime sweep!(state)

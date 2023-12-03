#!/usr/bin/julia

include("util.jl")

using BenchmarkTools

const param = Params{Float64}(1.0, 0.01, 1.0, 110, (24, 24, 24))
const state = State(param)

@time for i in 1:10000
    sweep!(state)
end

@time for i in 1:10000
    sweep!(state)
end

@btime for i in 1:1000
    sweep!(state)
end
@btime for i in 1:1000
    sweep!(state)
end

@time for i in 1:10000
    sweep!(state)
end

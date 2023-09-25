#!/usr/bin/julia

include("process_data.jl")

for f in readdir(joinpath(@__DIR__, "data"), join=true)
    bin = read_bin_compressed(f)
    @show f, count_cycles(bin)
end

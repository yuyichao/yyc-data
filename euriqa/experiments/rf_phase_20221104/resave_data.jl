#!/usr/bin/julia

using DelimitedFiles
using HDF5

function load_data(fname, tscale, vscale)
    data = readdlm(fname, ',', Float64, skipstart=5)
    tstart = round.(Int32, data[1, 1] .* tscale)
    tend = round.(Int32, data[end, 1] .* tscale)
    @assert size(data, 1) == (tend - tstart + 1)
    return Dict("time"=>[tstart:tend;],
                "value"=>round.(Int16, data[:, 2] .* vscale),
                "tscale"=>tscale, "vscale"=>vscale)
end

function resave_data(fname, data)
    h5open(fname, "w") do io
        for (k, v) in data
            io[k, compress=9] = v
        end
    end
end

resave_data(ARGS[2], load_data(ARGS[1], 2e9, 300.0))

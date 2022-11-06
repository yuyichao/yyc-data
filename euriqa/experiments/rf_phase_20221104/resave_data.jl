#!/usr/bin/julia

using DelimitedFiles
using HDF5

function guess_tscale(ts)
    dt = (ts[end] - ts[1]) / (length(ts) - 1)
    if 0.99e-9 <= dt <= 1.01e-9
        return 1e9
    elseif 0.49e-9 <= dt <= 0.51e-9
        return 2e9
    elseif 0.19e-9 <= dt <= 0.21e-9
        return 5e9
    end
    @show dt
    error()
end

function load_data(fname, vscale)
    data = readdlm(fname, ',', Float64, skipstart=5)
    tscale = guess_tscale(data[:, 1])
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

resave_data(ARGS[2], load_data(ARGS[1], 300.0))

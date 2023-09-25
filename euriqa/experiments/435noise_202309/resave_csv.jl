#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData.Agilent

function resave_csv(inname, outname, nwaveforms=nothing)
    waveforms = Agilent.load_bin(inname)
    if nwaveforms === nothing
        nwaveforms = length(waveforms)
    else
        @assert length(waveforms) >= nwaveforms
    end
    dt = round(waveforms[1].meta["XIncrement"]::Float64, digits=12)
    npoints = length(waveforms[1].data[1])
    for i in 1:nwaveforms
        @assert length(waveforms[i].data) == 1
        @assert eltype(waveforms[i].data[1]) === Float32
        @assert length(waveforms[i].data[1]) === npoints
    end
    midpoint = (npoints ÷ 2)
    times = ((-midpoint):(npoints - midpoint - 1)) .* dt

    open(outname, "w") do io
        for i in 1:npoints
            print(io, "$(times[i])")
            for j in 1:nwaveforms
                print(io, ",$(waveforms[j].data[1][i])")
            end
            println(io)
        end
    end
end

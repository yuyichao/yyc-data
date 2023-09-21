#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData.Agilent

function resave_data(inname, outname)
    waveforms = Agilent.load_bin(inname)
    @assert length(waveforms) == 3
    @assert length(waveforms[1].data) == 1
    @assert length(waveforms[2].data) == 1
    @assert length(waveforms[3].data) == 1

    @assert eltype(waveforms[1].data[1]) === Float32
    @assert eltype(waveforms[2].data[1]) === Float32

    vpoints = sort!(union(waveforms[1].data[1], waveforms[2].data[1]))
    vdeltas = sort!(unique!(diff(vpoints)))

    vdelta_threshold = max(vdeltas[end] * 1e-3, 0.0005)

    ifirst = searchsortedfirst(vdeltas, vdelta_threshold)

    delta1 = vdeltas[ifirst]

    ilast = searchsortedlast(vdeltas, delta1 * 1.2)

    delta2 = vdeltas[ilast]

    delta_avg = (delta1 + delta2) / 2

    cur_int_val = 0
    cur_float_val = vpoints[1]

    point_map = Dict{Float32,Int}()

    for v in vpoints
        d = v - cur_float_val
        if d <= 0.01 * delta_avg
            point_map[v] = cur_int_val
            continue
        end
        df, di = modf(d / delta_avg)
        di = round(Int, di)
        if df >= 0.8
            di += 1
        elseif df > 0.2
            error("Value spacing ($(d)) is not an integer multiple of step size ($(delta_avg)).")
        end
        cur_float_val = v
        cur_int_val += di
        point_map[v] = cur_int_val
    end

    int_offset = (cur_int_val + 1) ÷ 2

    if cur_int_val < 2^8
        ET = Int8
    elseif cur_int_val < 2^16
        ET = Int16
    elseif cur_int_val < 2^32
        ET = Int32
    else
        ET = Int64
    end

    dy = (vpoints[end] - vpoints[1]) / cur_int_val
    y0 = vpoints[1] + dy * int_offset

    for i in 1:2
        open("$(outname)_$(i).bin", "w") do fd
            waveform = waveforms[i]
            write(fd, Float64[waveform.meta["XOrigin"]::Float64,
                              waveform.meta["XIncrement"]::Float64, y0, dy])
            raw_data = waveform.data[1]
            write(fd, length(raw_data) % Int64)
            write(fd, sizeof(ET) % Int64)
            write(fd, ET[point_map[v] - int_offset for v in raw_data])
        end
    end
end

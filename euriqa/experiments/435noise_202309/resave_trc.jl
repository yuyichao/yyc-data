#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData.Lecroy

function resave_data(inname, outname)
    data = Lecroy.load(inname)
    @assert data.y2 === nothing
    x0 = data.x0[]
    dx = Float64(data.WAVEDESC["HORIZ_INTERVAL"]::Float32)
    y0 = -Float64(data.WAVEDESC["VERTICAL_OFFSET"]::Float32)
    dy = -Float64(data.WAVEDESC["VERTICAL_GAIN"]::Float32)
    y1 = data.y1
    @assert ndims(y1) == 1
    if eltype(y1) == Int16
        if all(y->(y % Int8 == 0), y1)
            y1 = Int8.(y1 .>> 8)
            dy *= 256
        end
    end
    coord = [x0, dx, y0, dy]
    open(outname, "w") do fd
        write(fd, coord)
        write(fd, length(y1) % Int64)
        write(fd, sizeof(eltype(y1)) % Int64)
        write(fd, y1)
    end
end

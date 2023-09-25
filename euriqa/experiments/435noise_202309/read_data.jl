#!/usr/bin/julia

using LibArchive

function read_bin_compressed(fname)
    LibArchive.Reader(fname) do reader
        LibArchive.support_format_raw(reader)
        LibArchive.support_filter_all(reader)
        LibArchive.next_header(reader)
        x0 = read(reader, Float64)
        dx = read(reader, Float64)
        y0 = read(reader, Float64)
        dy = read(reader, Float64)
        len = read(reader, Int64)
        elsz = read(reader, Int64)
        if elsz == 1
            res = read!(reader, Vector{Int8}(undef, len))
        else
            @assert elsz == 2
            res = read!(reader, Vector{Int16}(undef, len))
        end
        return (x0=x0, dx=dx, y0=y0, dy=dy, data=res)
    end
end

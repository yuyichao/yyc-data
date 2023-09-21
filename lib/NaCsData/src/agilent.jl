#

module Agilent

@inline function read_header(io, ::Type{T}, bytes_left) where T
    res = read(io, T)
    bytes_left[] -= sizeof(T)
    return res
end

@inline function read_header_str(io, count, bytes_left)
    res = String(read(io, count))
    bytes_left[] -= count
    return res
end

const _DataVectorType = Union{Vector{Float32},Vector{Int32},Vector{UInt8}}

struct Waveform
    meta::Dict{String,Any}
    times::typeof(1.0:0.1:1.1)
    data::Vector{_DataVectorType}
    Waveform(meta, times, data) = new(meta, times, data)
end

function load_bin(io::IO)
    # read file header
    fileCookie = read(io, 2)
    fileVersion = read(io, 2)
    fileSize = Int(read(io, Int32))
    nWaveforms = Int(read(io, Int32))

    # verify cookie
    if fileCookie != b"AG"
        error("Unrecognized file format.")
    end

    waveforms = Waveform[]

    for waveformIndex in 1:nWaveforms
        # read waveform header
        headerSize = Int(read(io, Int32))
        bytesLeft = Ref(headerSize - 4)

        meta = Dict{String,Any}()

        meta["WaveformType"] = Int(read_header(io, Int32, bytesLeft))
        meta["NWaveformBuffers"] =
            NWaveformBuffers = Int(read_header(io, Int32, bytesLeft))
        meta["NPoints"] = NPoints = Int(read_header(io, Int32, bytesLeft))
        meta["Count"] = Int(read_header(io, Int32, bytesLeft))
        meta["XDisplayRange"] = read_header(io, Float32, bytesLeft)
        meta["XDisplayOrigin"] = read_header(io, Float64, bytesLeft)
        meta["XIncrement"] = XIncrement = read_header(io, Float64, bytesLeft)
        meta["XOrigin"] = XOrigin = read_header(io, Float64, bytesLeft)
        meta["XUnits"] = Int(read_header(io, Int32, bytesLeft))
        meta["YUnits"] = Int(read_header(io, Int32, bytesLeft))
        meta["DateString"] = read_header_str(io, 16, bytesLeft)
        meta["TimeString"] = read_header_str(io, 16, bytesLeft)
        meta["FrameString"] = read_header_str(io, 24, bytesLeft)
        meta["WaveformString"] = read_header_str(io, 16, bytesLeft)
        meta["TimeTag"] = read_header(io, Float64, bytesLeft)
        meta["SegmentIndex"] = Int(read_header(io, UInt32, bytesLeft))
        # skip over any remaining data in the header
        if bytesLeft[] < 0
            error("Header size smaller than expected.")
        end
        skip(io, bytesLeft[])
        # generate time vector from xIncrement and xOrigin values
        times = XIncrement .* (0:NPoints - 1) .+ XOrigin
        waveform = Waveform(meta, times,
                            Vector{_DataVectorType}(undef, NWaveformBuffers))
        push!(waveforms, waveform)

        for bufferIndex in 1:NWaveformBuffers
            # read waveform buffer header
            headerSize = Int(read(io, Int32))
            bytesLeft = Ref(headerSize - 4)

            BufferType = read_header(io, Int16, bytesLeft)
            BytesPerPoint = read_header(io, Int16, bytesLeft)
            BufferSize = read_header(io, Int32, bytesLeft)

            # skip over any remaining data in the header
            if bytesLeft[] < 0
                error("Header size smaller than expected.")
            end
            skip(io, bytesLeft[])

            if ((BufferType == 1) | (BufferType == 2) | (BufferType == 3))
                # bufferType is PB_DATA_NORMAL, PB_DATA_MIN, or PB_DATA_MAX (float)
                waveform.data[bufferIndex] =
                    read!(io, Vector{Float32}(undef, NPoints))
            elseif (bufferType == 4)
                # bufferType is PB_DATA_COUNTS (int32)
                waveform.data[bufferIndex] =
                    read!(io, Vector{Int32}(undef, NPoints))
            elseif (bufferType == 5)
                # bufferType is PB_DATA_LOGIC (int8)
                waveform.data[bufferIndex] =
                    read!(io, Vector{UInt8}(undef, NPoints))
            else
                # unrecognized bufferType read as unformated bytes
                waveform.data[bufferIndex] =
                    read!(io, Vector{UInt8}(undef, BufferSize))
            end
        end
    end
    return waveforms
end

load_bin(filename::AbstractString) = open(load_bin, filename)

end

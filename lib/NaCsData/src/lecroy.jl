#

module Lecroy

# This is a module for reading binary files from LeCroy scopes.
#
# LeCroy binary files have a .trc extension. This module reads version 'LECROY_2_3' of
# the format (an exception is raised if a different format is encountered).
#
# Based on lecroy.py
#
# D. Guarisco, 2013-2018. Assembled from various sources.
#
# Version history:
#   1.0         2013-02-10  First release
#   1.1         2013-02-13  Added support for sequence acquisitions. Fixed a few bugs.
#   1.2         2013-07-30  Correctly handles the case where WAVEDESC block starts at the                            beginning of the file (i.e., there is no file size info)
#   1.3         2013-08-27  Correctly imports FFT traces. TIMEBASE and FIXED VERT GAIN
#                           reflect horizontal and vertical units. 
#   1.4         2018-12-24  Ported to Python 3 (no longer works in Python 2). 
#                           Correctly reports FILE_SIZE in WAVEDESC.
#   1.5         2019-04-27  Replaced string decoding from standard 'utf-8' to 'latin_1'
#                           because utf-8 would generate errors on some .trc files. 
#                           Fixed bug with file size (try to read it from file and rely on
#                           OS if it's missing).
#

@enum Endian LittleEndian BigEndian

@enum RecordType begin
    SingleSweep = 0
    Interleaved = 1
    Histogram = 2
    Graph = 3
    FilterCoefficient = 4
    Complex = 5
    Extrema = 6
    SequenceObsolete = 7
    CenteredRIS = 8
    PeakDetect = 9
end

@enum ProcessingDone begin
    NoProcessing = 0
    FIRFilter = 1
    Interpolated = 2
    Sparsed = 3
    Autoscaled = 4
    NoResult = 5
    Rolling = 6
    Cumulative = 7
end

@enum VertCoupling begin
    DC50Ohms = 0
    Ground = 1
    DC1MOhm = 2
    # Ground = 3
    AC1MOhm = 4
end

function read_fixed_len_str(io, len)
    buf = IOBuffer()
    ended = false
    for i in 1:len
        c = read(io, UInt8)
        if c == 0
            ended = true
        end
        if !ended
            write(buf, c)
        end
    end
    return String(take!(buf))
end

"""
Reads a binary LeCroy file (file extension: .trc). Only version 'LECROY_2_3' of the
format is supported. A detailed description of the format can be obtained from the
instrument by issuing the command 'TEMPLATE?'

The following results are returned:
    (WAVEDESC,USER_TEXT,x,y1,y2), where

WAVEDESC is the file header, containing information about the waveform
USER_TEXT is an optional string (empty by default)
x  is a numpy array (dtype='float64') containing the time offsets of the waveform
   relative to the trigger time
y1 is a numpy array (dtype='float64') containing the primary waveform
y2 is a numpy array (dtype='float64') containing the secondary waveform (empty for a
   single sweep acquisition).

The shape of the arrays x and y depends on the mode of the scope during the data
acquisition. If the scope is not in sequence mode (SUBARRAY_COUNT=1) a single sweep
was captured. In this case, x and y1 are one-dimensional arrays of size
WAVE_ARRAY_COUNT. If the scope was in sequence mode, several sweeps (SUBARRAY_COUNT)
were captured in sequence. In this case, the shape of x and y is a two-dimensional
array(SUBARRAY_COUNT,WAVE_ARRAY_COUNT/SUBARRAY_COUNT).

    Parameter name                          Type
    DESCRIPTOR_NAME                         String
    TEMPLATE_NAME                           String
    COMM_TYPE                               Endian
    COMM_ORDER                              Int (0 - byte, 1 - word)
    WAVE_DESCRIPTOR                         Int32
    USER_TEXT                               Int32
    RES_DESC1                               Int32
    TRIGTIME_ARRAY                          Int32
    RIS_TIME_ARRAY                          Int32
    RES_ARRAY1                              Int32
    WAVE_ARRAY_1                            Int32
    WAVE_ARRAY_2                            Int32
    RES_ARRAY2                              Int32
    RES_ARRAY3                              Int32
    INSTRUMENT_NAME                         String
    INSTRUMENT_NUMBER                       Int32
    TRACE_LABEL                             String
    RESERVED1                               Int16
    RESERVED2                               Int16
    WAVE_ARRAY_COUNT                        Int32
    PNTS_PER_SCREEN                         Int32
    FIRST_VALID_PNT                         Int32
    LAST_VALID_PNT                          Int32
    FIRST_POINT                             Int32
    SPARSING_FACTOR                         Int32
    SEGMENT_INDEX                           Int32
    SUBARRAY_COUNT                          Int32
    SWEEPS_PER_ACQ                          Int32
    POINTS_PER_PAIR                         Int16
    PAIR_OFFSET                             Int16
    VERTICAL_GAIN                           Float32
    VERTICAL_OFFSET                         Float32
    MAX_VALUE                               Float32
    MIN_VALUE                               Float32
    NOMINAL_BITS                            Int16
    NOM_SUBARRAY_COUNT                      Int16
    HORIZ_INTERVAL                          Float32
    HORIZ_OFFSET                            Float64
    PIXEL_OFFSET                            Float64
    VERTUNIT                                String
    HORUNIT                                 String
    HORIZ_UNCERTAINTY                       Float32
    ACQ_DURATION                            Float32
    RECORD_TYPE                             RecordType
    PROCESSING_DONE                         ProcessingDone
    RESERVED5                               Int16
    RIS_SWEEPS                              Int16
    TIMEBASE                                Int (0 for external)
    VERT_COUPLING                           VertCoupling
    PROBE_ATT                               Float32
    FIXED_VERT_GAIN                         Int
    BANDWIDTH_LIMIT                         Bool
    VERTICAL_VERNIER                        Float32
    ACQ_VERT_OFFSET                         Float32
    WAVE_SOURCE                             Int (1-based)
    FILE_SIZE                               Union{Int,Nothing}
"""
function load(io)
    read(io, UInt8) # "#"
    read(io, UInt8) # "9": number of digits representing the file size
    file_size_str = readuntil(io, "WAVEDESC")
    file_size = tryparse(Int, file_size_str)
    DESCRIPTOR_NAME = "WAVEDESC" * readuntil(io, '\0')
    for i in 1:(15 - length(DESCRIPTOR_NAME))
        read(io, UInt8) # "\0"
    end
    TEMPLATE_NAME = read_fixed_len_str(io, 16)
    comm_type1 = read(io, UInt8)
    comm_type2 = read(io, UInt8)

    comm_order1 = read(io, UInt8)
    comm_order2 = read(io, UInt8)
    if comm_order1 == comm_order2 == 0x0
        endian = BigEndian
        comm_type = (UInt16(comm_type1) << 8) | comm_type2
    else
        endian = LittleEndian
        comm_type = (UInt16(comm_type2) << 8) | comm_type1
    end

    if comm_type != 0 && comm_type != 1
        error("Unknown data type: $(comm_type)")
    end

    WAVEDESC = Dict{String,Any}("DESCRIPTOR_NAME"=>DESCRIPTOR_NAME,
                                "TEMPLATE_NAME"=>TEMPLATE_NAME,
                                "COMM_ORDER"=>endian,
                                "COMM_TYPE"=>comm_type,
                                "FILE_SIZE"=>file_size)

    if endian == BigEndian
        return _load(io, Val(BigEndian), comm_type, WAVEDESC)
    else
        return _load(io, Val(LittleEndian), comm_type, WAVEDESC)
    end
end

function _read(io, ::Type{T}, ::Val{E}) where {T, E}
    v = read(io, T)
    if E == BigEndian
        return ntoh(v)
    else
        return ltoh(v)
    end
end

function _load(io, E, comm_type, WAVEDESC)
    WAVEDESC["WAVE_DESCRIPTOR"] = _read(io, Int32, E)
    WAVEDESC["USER_TEXT"] = USER_TEXT = _read(io, Int32, E)
    WAVEDESC["RES_DESC1"] = _read(io, Int32, E)
    WAVEDESC["TRIGTIME_ARRAY"] = _read(io, Int32, E)
    WAVEDESC["RIS_TIME_ARRAY"] = _read(io, Int32, E)
    WAVEDESC["RES_ARRAY1"] = _read(io, Int32, E)
    WAVEDESC["WAVE_ARRAY_1"] = _read(io, Int32, E)
    WAVEDESC["WAVE_ARRAY_2"] = WAVE_ARRAY_2 = _read(io, Int32, E)
    WAVEDESC["RES_ARRAY2"] = _read(io, Int32, E)
    WAVEDESC["RES_ARRAY3"] = _read(io, Int32, E)
    WAVEDESC["INSTRUMENT_NAME"] = read_fixed_len_str(io, 16)
    WAVEDESC["INSTRUMENT_NUMBER"] = _read(io, Int32, E)
    WAVEDESC["TRACE_LABEL"] = read_fixed_len_str(io, 16)
    WAVEDESC["RESERVED1"] = _read(io, Int16, E)
    WAVEDESC["RESERVED2"] = _read(io, Int16, E)
    WAVEDESC["WAVE_ARRAY_COUNT"] = WAVE_ARRAY_COUNT = _read(io, Int32, E)
    WAVEDESC["PNTS_PER_SCREEN"] = _read(io, Int32, E)
    WAVEDESC["FIRST_VALID_PNT"] = _read(io, Int32, E)
    WAVEDESC["LAST_VALID_PNT"] = _read(io, Int32, E)
    WAVEDESC["FIRST_POINT"] = _read(io, Int32, E)
    WAVEDESC["SPARSING_FACTOR"] = _read(io, Int32, E)
    WAVEDESC["SEGMENT_INDEX"] = _read(io, Int32, E)
    WAVEDESC["SUBARRAY_COUNT"] = SUBARRAY_COUNT = _read(io, Int32, E)
    WAVEDESC["SWEEPS_PER_ACQ"] = _read(io, Int32, E)
    WAVEDESC["POINTS_PER_PAIR"] = _read(io, Int16, E)
    WAVEDESC["PAIR_OFFSET"] = _read(io, Int16, E)
    WAVEDESC["VERTICAL_GAIN"] = _read(io, Float32, E)
    WAVEDESC["VERTICAL_OFFSET"] = _read(io, Float32, E)
    WAVEDESC["MAX_VALUE"] = _read(io, Float32, E)
    WAVEDESC["MIN_VALUE"] = _read(io, Float32, E)
    WAVEDESC["NOMINAL_BITS"] = _read(io, Int16, E)
    WAVEDESC["NOM_SUBARRAY_COUNT"] = _read(io, Int16, E)
    WAVEDESC["HORIZ_INTERVAL"] = _read(io, Float32, E)
    WAVEDESC["HORIZ_OFFSET"] = HORIZ_OFFSET = _read(io, Float64, E)
    WAVEDESC["PIXEL_OFFSET"] = _read(io, Float64, E)
    WAVEDESC["VERTUNIT"] = read_fixed_len_str(io, 48)
    WAVEDESC["HORUNIT"] = read_fixed_len_str(io, 48)
    WAVEDESC["HORIZ_UNCERTAINTY"] = _read(io, Float32, E)
    WAVEDESC["TRIGGER_TIME_SECONDS"] = _read(io, Float64, E)
    WAVEDESC["TRIGGER_TIME_MINUTES"] = read(io, UInt8)
    WAVEDESC["TRIGGER_TIME_HOURS"] = read(io, UInt8)
    WAVEDESC["TRIGGER_TIME_DAYS"] = read(io, UInt8)
    WAVEDESC["TRIGGER_TIME_MONTHS"] = read(io, UInt8)
    WAVEDESC["TRIGGER_TIME_YEAR"] = _read(io, Int16, E)
    WAVEDESC["TRIGGER_TIME_UNUSED"] = _read(io, Int16, E)
    WAVEDESC["ACQ_DURATION"] = _read(io, Float32, E)
    WAVEDESC["RECORD_TYPE"] = RecordType(_read(io, Int16, E))
    WAVEDESC["PROCESSING_DONE"] = ProcessingDone(_read(io, Int16, E))
    WAVEDESC["RESERVED5"] = _read(io, Int16, E)
    WAVEDESC["RIS_SWEEPS"] = _read(io, Int16, E)
    _TIMEBASE = _read(io, Int16, E)
    if _TIMEBASE == 100
        WAVEDESC["TIMEBASE"] = 0
    else
        mant = (1, 2, 5)[_TIMEBASE % 3 + 1]
        WAVEDESC["TIMEBASE"] = mant * exp10(_TIMEBASE ÷ 3 - 12)
    end
    _VERT_COUPLING = _read(io, Int16, E)
    if _VERT_COUPLING == 3
        _VERT_COUPLING = 1
    end
    WAVEDESC["VERT_COUPLING"] = VertCoupling(_VERT_COUPLING)
    WAVEDESC["PROBE_ATT"] = _read(io, Float32, E)
    _FIXED_VERT_GAIN = _read(io, Int16, E)
    mant = (1, 2, 5)[_FIXED_VERT_GAIN % 3 + 1]
    WAVEDESC["FIXED_VERT"] = mant * exp10(_FIXED_VERT_GAIN ÷ 3 - 6)
    WAVEDESC["BANDWIDTH_LIMIT"] = Bool(_read(io, Int16, E))
    WAVEDESC["VERTICAL_VERNIER"] = _read(io, Float32, E)
    WAVEDESC["ACQ_VERT_OFFSET"] = _read(io, Float32, E)
    WAVEDESC["WAVE_SOURCE"] = _read(io, Int16, E) + 1

    # Read user text (160 char. maximum)
    WAVEDESC["TEXT"] = String(read(io, USER_TEXT))

    # Read waveforms. Distinguish case of acquisition sequence or single acquisition
    y2 = nothing
    if SUBARRAY_COUNT > 1
        # Multiple segments
        # Sanity check
        if WAVE_ARRAY_COUNT % SUBARRAY_COUNT != 0
            error("Number of data points is not a multiple of number of segments")
        end
        npts = WAVE_ARRAY_COUNT ÷ SUBARRAY_COUNT
        # Read TRIGTIME array first.
        # There are SUBARRAY_COUNT repetitions of two doubles,
        # TRIGGER_TIME and TRIGGER_OFFSET
        x0 = Vector{Float64}(undef, SUBARRAY_COUNT)
        for i in 1:SUBARRAY_COUNT
            trig_time = _read(io, Float64, E)
            trig_offset = _read(io, Float64, E)
            x0[i] = trig_time + trig_offset
        end
        # Now read data array
        if comm_type == 0
            y1 = Matrix{Int8}(undef, npts, SUBARRAY_COUNT)
            read!(io, y1)
        else
            y1 = Matrix{Int16}(undef, npts, SUBARRAY_COUNT)
            read!(io, y1)
            if E == BigEndian
                y1 .= ntoh(y1)
            else
                y1 .= ltoh(y1)
            end
        end
    else
        # Single sweep. Read waveforms from file
        if comm_type == 0
            y1 = Vector{Int8}(undef, WAVE_ARRAY_COUNT)
            read!(io, y1)
            if WAVE_ARRAY_2 > 0
                y2 = Vector{Int8}(undef, WAVE_ARRAY_COUNT)
                read!(io, y2)
            end
        else
            y1 = Vector{Int16}(undef, WAVE_ARRAY_COUNT)
            read!(io, y1)
            if E == BigEndian
                y1 .= ntoh(y1)
            else
                y1 .= ltoh(y1)
            end
            if WAVE_ARRAY_2 > 0
                y2 = Vector{Int16}(undef, WAVE_ARRAY_COUNT)
                read!(io, y2)
                if E == BigEndian
                    y2 .= ntoh(y2)
                else
                    y2 .= ltoh(y2)
                end
            end
            # Generate time intervals
            x0 = [HORIZ_OFFSET]
        end
    end
    return WAVEDESC, x0, y1, y2
end

end

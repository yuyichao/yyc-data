#

include("am_shared.jl")
include("ind-modulation.jl")

using GoldGates
import ProtoBuf as PB

params_file = "072125_goldparams_13ions.json"
sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
end

ωs = 2π .* sysparams.modes.radial1

candidates_file = "data_am/am_candidates_20260522.binpb"
candidates = open(candidates_file) do io
    decoder = PB.ProtoDecoder(io)
    candidates = PB.decode(decoder, Candidates)
    return candidates.candidates
end

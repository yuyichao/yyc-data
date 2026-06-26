#

include("am_shared.jl")

using GoldGates
using LinearAlgebra

import ProtoBuf as PB

params_file = "072125_goldparams_13ions.json"

sysparams = open(params_file) do io
    read(io, GoldGates.SystemParams; format=:json)
end
ωs = 2π .* sysparams.modes.radial1
nions = length(ωs)
modes = Seq.Modes()
for ω in ωs
    push!(modes, ω)
end

candidates_file = ARGS[1]
candidates = open(candidates_file) do io
    decoder = PB.ProtoDecoder(io)
    candidates = PB.decode(decoder, Candidates)
    return candidates.candidates
end

ncands = length(candidates)
NSeg = length(candidates[1].Ωs) + 1
τ = candidates[1].τ
ω = candidates[1].ω

cand_arys = [candidates[i].Ωs[j] for j in 1:NSeg - 1, i in 1:ncands]

F = svd(cand_arys)

cand_ary_filter = F.U[:, F.S .>= F.S[1] * 0.01]

freq_spec = Seq.FreqSpec(false, sym=true)
amp_spec = Seq.AmpSpec(cb=get_am_cbs(NSeg), sym=false)

pre_obj = avg_area_obj(NSeg, modes, SL.pmask_full, freq=freq_spec, amp=amp_spec)
s = Seq.Summarizer{NSeg}()

candidates_filtered = Candidate[]

for i in 1:size(cand_ary_filter, 2)
    args = [τ; cand_ary_filter[:, i]; ω]
    raw_params = Seq.RawParams(pre_obj, args)
    props = get(s, raw_params, modes)
    push!(candidates_filtered, Candidate(args, pre_obj.param, props))
end

open(candidates_file, "w") do io
    encoder = PB.ProtoEncoder(io)
    PB.encode(encoder, Candidates(candidates_filtered, ""))
end

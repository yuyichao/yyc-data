#

import MSSim: Optimizers as Opts, SegSeq as SS, SymLinear as SL, Sequence as Seq, Utils as U
using StaticArrays

function get_am_cbs(NSeg)
    return ntuple(i -> begin
                      prev = (i - 1) / (NSeg / 2) - 1
                      mid = i / (NSeg / 2) - 1
                      next = (i + 1) / (NSeg / 2) - 1
                      function (x)
                          if x < prev || x > next
                              return 0.0
                          elseif x < mid
                              return (x - prev) / (mid - prev)
                          else
                              return (next - x) / (next - mid)
                          end
                      end
                  end, NSeg - 1)
end

struct AvgAreaObjCallback{NModes}
    dis_weight::Float64
    disδ_weight::Float64
    area_weights::MVector{NModes,Float64}
    function AvgAreaObjCallback(NModes, dis_weight, disδ_weight, area_weights)
        cb = new{NModes}(dis_weight, disδ_weight, MVector{NModes,Float64}(undef))
        cb.area_weights .= area_weights
        return cb
    end
end

function (obj::AvgAreaObjCallback{NModes})(vals) where NModes
    @assert length(vals) == NModes * 3 + 1
    weights = obj.area_weights
    @inbounds begin
        v1 = muladd(vals[1], obj.dis_weight, vals[1 + NModes] * obj.disδ_weight)
        v2 = abs(vals[1 + NModes * 2]) * weights[1]
    end
    @inbounds @simd ivdep for i in 2:NModes
        v1 = muladd(vals[i], obj.dis_weight, muladd(vals[i + NModes], obj.disδ_weight, v1))
        v2 = muladd(abs(vals[i + NModes * 2]), weights[i], v2)
    end
    v1 = v1 + 1e-5
    return v1 / v2 * vals[end]
end

function avg_area_obj(nseg, modes, pmask;
                      freq=Seq.FreqSpec(), amp=Seq.AmpSpec(),
                      dis_weight=5, disδ_weight=0.01)
    nmodes = length(modes.modes)
    dis_args = ntuple(i->(:dis2, i), nmodes)
    disδ_args = ntuple(i->(:disδ2, i), nmodes)
    area_args = ntuple(i->(:area, i), nmodes)

    mask_dis_area = SS.ValueMask(true, true, true, false, true, false)
    buf = SL.ComputeBuffer{nseg,Float64}(Val(SS.mask_allδ), Val(SS.mask_allδ))

    area_weights = zeros(nmodes)
    area_weights[1] = 1
    area_weights[2] = 1
    return Seq.Objective(pmask, (dis_args..., disδ_args..., area_args..., (:τ, 0)),
                         Opts.autodiff(AvgAreaObjCallback(nmodes, dis_weight,
                                                          disδ_weight, area_weights)),
                         modes, buf, freq=freq, amp=amp)
end

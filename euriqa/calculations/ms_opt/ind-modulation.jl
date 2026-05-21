#!/usr/bin/julia

using MSSim: Utils as U

const rf_vv = U.cos_f1
const rf_vd = U.sin_f1
const rf_dv = U.TrigRatio{true,3,(),(1,),(-1,)}()
const rf_dd = U.cos_f3

const if_vv = U.TrigRatio{true,2,(1,),(-1,),()}()
const if_vd = U.TrigRatio{false,3,(-1,1//2),(),(1,)}()
const if_dv = U.TrigRatio{false,3,(1,1//2),(-1,),(-1,)}()
const if_dd = U.sin_f3

@inline function enclosed_area2_kernel(o1, o1′, o2, o2′, d, s, c)
    c_vv = @inline rf_vv(d, s, c)
    c_vd = @inline rf_vd(d, s, c)
    c_dv = @inline rf_dv(d, s, c)
    c_dd = @inline rf_dd(d, s, c)

    vv = o1 * o2
    vd = o1 * o2′
    dv = o1′ * o2
    dd = o1′ * o2′

    return muladd(c_vv,  vv, muladd(c_vd, vd, muladd(c_dv, dv, c_dd * dd)))
end

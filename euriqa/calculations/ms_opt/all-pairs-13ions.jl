#!/usr/bin/julia

include("opt-all-pair.jl")

nseg = 200
η_Yb171(f) = 0.1924385393508112 / sqrt(f)
const fs = [2.29049, 2.31223, 2.33447, 2.35660, 2.37964, 2.40187, 2.42345,
            2.44437, 2.46325, 2.48068, 2.49643, 2.50981, 2.51649]
const ηs = η_Yb171.(fs)
const bij = [0.002163619698947626 -0.019894333993010195 0.07990263247850947 -0.1959542649546874 0.3409659166835623 -0.45954763589500497 0.5047281319628083 -0.4595476358947884 0.3409659166837989 -0.19595426495458965 0.07990263247851344 -0.019894333993005258 0.0021636196989478743
             0.0071396181636678106 -0.05753283102846898 0.19547357544098015 -0.37938832365502567 0.46052353394342116 -0.3200182105272044 0.0 0.3200182105273888 -0.46052353394369755 0.3793883236548021 -0.1954735754409945 0.05753283102846495 -0.007139618163669534
             -0.018249239558250613 0.12361462840018689 -0.33254420091786313 0.44472648944450016 -0.22872005803180523 -0.19554847300587097 0.4134417073373916 -0.19554847300472525 -0.22872005803156883 0.4447264894437476 -0.33254420091756615 0.1236146284000512 -0.01824923955823083
             -0.0398846241690199 0.21868897045586588 -0.43028589455025196 0.29215796493448437 0.1863175039870539 -0.38127144023462645 0.0 0.38127144023514115 -0.18631750398830155 -0.2921579649354758 0.4302858945521276 -0.21868897045675953 0.039884624169179705
             0.07712220007227492 -0.32780345715192233 0.41319645824438495 0.032076667328037364 -0.37977444162142493 -0.002701578689309821 0.37576830363562647 -0.002701578689492867 -0.3797744416199687 0.03207666732788457 0.4131964582420578 -0.3278034571499661 0.07712220007182097
             -0.13409522479694663 0.4177755684963474 -0.2434744870926441 -0.31655881037945377 0.144563085647049 0.3565028092834569 0.0 -0.35650280928392075 -0.14456308564819773 0.316558810380403 0.24347448709360936 -0.417775568497799 0.13409522479741495
             -0.21125301253887932 0.44419275395368024 0.03180655196179913 -0.33738373307931335 -0.24588559566306237 0.1413373146118317 0.35437144150833166 0.14133731461085156 -0.245885595662734 -0.3373837330788401 0.03180655196168993 0.4441927539533746 -0.2112530125387379
             -0.30217912088262994 0.37043411370382284 0.2814435836968718 -0.07002662849707697 -0.32491388364205287 -0.2859774379576588 0.0 0.2859774379578876 0.32491388364227425 0.07002662849688487 -0.2814435836968081 -0.37043411370377427 0.302179120882572
             -0.3914400374056656 0.19247382191126752 0.35988979019519046 0.2517087393320628 0.00013533283088761348 -0.2421844792450492 -0.3411663352374448 -0.24218447924483172 0.00013533283060922018 0.25170873933215426 0.3598897901952505 0.19247382191128018 -0.3914400374056827
             -0.4554774750079228 -0.044926026877461854 0.211605450619302 0.3303919021675469 0.3164620569467598 0.19088743445490436 0.0 -0.19088743445506165 -0.3164620569466799 -0.3303919021674821 -0.2116054506193178 0.04492602687748547 0.45547747500789554
             -0.46790836191227264 -0.2561095095683251 -0.06959321731844297 0.0931837610333706 0.22319930737608462 0.3082480026043195 0.33796003557071747 0.30824800260425467 0.22319930737616317 0.09318376103329376 -0.06959321731841982 -0.25610950956833467 -0.4679083619122754
             0.42615110225953595 0.36754442103650975 0.30557641737122976 0.23756305532151642 0.16305861335609598 0.08309676955971396 0.0 -0.08309676955974268 -0.16305861335598906 -0.23756305532148866 -0.305576417371162 -0.36754442103645346 -0.4261511022594881
             0.2773500981126128 0.2773500981125969 0.27735009811258304 0.2773500981125812 0.2773500981125842 0.27735009811259226 0.277350098112603 0.277350098112613 0.277350098112625 0.2773500981126366 0.277350098112647 0.27735009811265454 0.27735009811265965]

const prefix = ARGS[1]
const nrep = parse(Int, ARGS[2])

const amp_ratio = 0.7
meta = Dict("amp_ratio"=>amp_ratio, "nseg"=>nseg)
_meta, candidates = load_candidates_dir(dirname(prefix))
println("Loaded $(length(candidates))")
@assert _meta === nothing || _meta == meta
const pre_pool = ThreadObjectPool() do
    return PreOptimizer{nseg}(2π .* fs, 2π .* (fs .+ 0.3); amp_ratio=amp_ratio,
                              tmin=250, tmax=400, ωmin=2π * 2.28, ωmax=2π * 2.497)
end
candidates = @time opt_all_rounds!(pre_pool, nrep, candidates)
@show length(candidates)
save_candidates(prefix, candidates, meta)

# full_opt = Optimizer{nseg}(2π .* fs, ηs, bij, 2π .* (fs .+ 0.3);
#                            tmin=250, tmax=400, ωmin=2π * 2.28, ωmax=2π * 2.3978)

# preopt_all_rounds!(full_opt, 10)
# @show length(full_opt.candidates)

# function (obj::ARobustObjective{NModes})(vals) where NModes
#     dis = vals[1]
#     disδ = vals[2]
#     area = zero(eltype(vals))
#     areaδ = zero(eltype(vals))
#     for i in 1:NModes
#         weight = obj.weights[i]
#         areai = vals[2 + i] * weight
#         areaδi = vals[2 + NModes + i] * weight
#         area += areai
#         areaδ += areaδi
#     end
#     return 5 * dis + disδ / 100 + (abs(area) - π / 2)^2 * 100 + (areaδ / 1e4)^2
# end

# struct ParamInfo
#     params::Vector{Float64}
#     dis::Float64
#     disδ::Float64
#     area::Float64
#     areaδ::Float64
#     total_t::Float64
#     maxΩ::Float64
# end

# function optimize_pairs(kernel, obj, preopt)
#     pairs = Set{Tuple{Int,Int}}()
#     for ion1 in 3:6
#         for ion2 in 2:ion1 - 1
#             push!(pairs, (ion1, ion2))
#         end
#     end

#     candidates_checked = 0
#     function check_modes(candidates_checked)
#         for (ion1, ion2) in collect(pairs)
#             mode_weight!(kernel.weights, ion1, ion2)
#             max_area = 0.0
#             found = false
#             for param_idx in candidates_checked + 1:length(candidates)
#                 params = candidates[param_idx]
#                 dis, disδ, area, areaδ, total_t, maxΩ = compute_properties(params)
#                 max_area = max(max_area, abs(area))
#                 if abs(area) >= 1.56
#                     println("Accept $(abs(area)) for $(ion1), $(ion2)")
#                     found = true
#                     delete!(pairs, (ion1, ion2))
#                     break
#                 end
#             end
#             if !found
#                 println("Reject $(max_area) for $(ion1), $(ion2)")
#             end
#         end
#         return length(candidates)
#     end
#     candidates = Vector{Float64}[]
#     @time pre_optimize!(preopt, 500, candidates)
#     while !isempty(pairs)
#         @time pre_optimize!(preopt, 100, candidates)
#         candidates_checked = check_modes(candidates_checked)
#         @show candidates_checked, pairs
#     end
#     pre_optimize!(preopt, 500, candidates)

#     nargs = Seq.nparams(obj)
#     tracker = Opts.NLVarTracker(nargs)
#     Opts.set_bound!(tracker, obj.param.τ, 0.1, 3)
#     Opts.set_bound!(tracker, obj.param.Ωs[1], 0.01, 0.41)

#     # opt = NLopt.Opt(:LD_SLSQP, nargs) # 13 ms
#     opt = NLopt.Opt(:LD_LBFGS, nargs) # 50 ms
#     # NLopt.xtol_rel!(opt, 1e-7)
#     # NLopt.ftol_rel!(opt, 1e-7)
#     NLopt.min_objective!(opt, obj)
#     NLopt.lower_bounds!(opt, Opts.lower_bounds(tracker))
#     NLopt.upper_bounds!(opt, Opts.upper_bounds(tracker))
#     NLopt.maxeval!(opt, 10000)
#     for ion1 in 3:6
#         for ion2 in 2:ion1 - 1
#             mode_weight!(kernel.weights, ion1, ion2)
#             tier1 = ParamInfo[]
#             tier2 = ParamInfo[]
#             tier3 = ParamInfo[]
#             tier4 = ParamInfo[]
#             for params in candidates
#                 dis, disδ, area, areaδ, total_t, maxΩ = compute_properties(params)
#                 if abs(area) < 1.56
#                     continue
#                 end
#                 info = ParamInfo(params, dis, disδ, area, areaδ, total_t, maxΩ)
#                 if abs(areaδ) < abs(area) * 2 && abs(area) > 1.8 && total_t <= 220
#                     push!(tier1, info)
#                 elseif abs(areaδ) < abs(area) * 2 && (abs(area) > 1.8 || total_t <= 220)
#                     push!(tier2, info)
#                 elseif abs(area) > 1.8 || total_t <= 220
#                     push!(tier3, info)
#                 else
#                     push!(tier4, info)
#                 end
#             end
#             if !isempty(tier1)
#                 mode_cand = tier1
#             elseif !isempty(tier2)
#                 mode_cand = tier2
#             elseif !isempty(tier3)
#                 mode_cand = tier3
#             else
#                 @assert !isempty(tier4)
#                 mode_cand = tier4
#             end
#             println("Found $(length(mode_cand)) candidates for $ion1, $ion2")
#             # sort!(mode_cand, by=x->(abs(x.areaδ / x.area), 1 / abs(x.area), x.total_t))

#             min_areaδ = 0.0
#             max_area = 0.0
#             min_total_t = 0.0
#             for info in mode_cand
#                 areaδ = abs(info.areaδ / info.area)
#                 area = abs(info.area)
#                 total_t = info.total_t
#                 min_areaδ = min(areaδ, min_areaδ)
#                 max_area = max(area, max_area)
#                 min_total_t = min(total_t, min_total_t)
#             end
#             tier1 = ParamInfo[]
#             tier2 = ParamInfo[]
#             tier3 = ParamInfo[]
#             tier4 = ParamInfo[]
#             for info in mode_cand
#                 areaδ = abs(info.areaδ / info.area)
#                 area = abs(info.area)
#                 total_t = info.total_t
#                 if (areaδ < min_areaδ * 3 &&
#                     (area > max_area * 0.7 && total_t <= min_total_t + 20))
#                     push!(tier1, info)
#                 elseif (areaδ < min_areaδ * 3 &&
#                     (area > max_area * 0.7 || total_t <= min_total_t + 20))
#                     push!(tier2, info)
#                 elseif area > max_area * 0.7 || total_t <= min_total_t + 20
#                     push!(tier3, info)
#                 else
#                     push!(tier4, info)
#                 end
#             end
#             if !isempty(tier1)
#                 mode_cand = tier1
#             elseif !isempty(tier2)
#                 mode_cand = tier2
#             elseif !isempty(tier3)
#                 mode_cand = tier3
#             else
#                 @assert !isempty(tier4)
#                 mode_cand = tier4
#             end
#             println("Found $(length(mode_cand)) candidates for $ion1, $ion2")

#             tier1 = ParamInfo[]
#             tier2 = ParamInfo[]
#             tier3 = ParamInfo[]
#             tier4 = ParamInfo[]
#             for old_info in mode_cand
#                 params = copy(old_info.params)
#                 params[obj.param.Ωbase] *= sqrt(π / 2 / abs(old_info.area))
#                 objv, params, ret = NLopt.optimize!(opt, params)
#                 dis, disδ, area, areaδ, total_t, maxΩ = compute_properties(params)
#                 info = ParamInfo(params, dis, disδ, area, areaδ, total_t, maxΩ)
#                 if abs(areaδ) < abs(area) * 1.5 && maxΩ < 0.6 && total_t <= 220
#                     push!(tier1, info)
#                 elseif abs(areaδ) < abs(area) * 1.5 && (maxΩ < 0.6 || total_t <= 220)
#                     push!(tier2, info)
#                 elseif maxΩ < 0.6 || total_t <= 220
#                     push!(tier3, info)
#                 else
#                     push!(tier4, info)
#                 end
#             end
#             if !isempty(tier1)
#                 mode_cand = tier1
#             elseif !isempty(tier2)
#                 mode_cand = tier2
#             elseif !isempty(tier3)
#                 mode_cand = tier3
#             else
#                 @assert !isempty(tier4)
#                 mode_cand = tier4
#             end
#             println("Found $(length(mode_cand)) candidates for $ion1, $ion2")
#             sort!(mode_cand, by=x->(abs(x.areaδ / x.area), 1 / abs(x.area), x.total_t))
#             info0 = mode_cand[1]
#             @show info0.dis, info0.disδ, info0.area, info0.areaδ, info0.total_t, info0.maxΩ
#         end
#     end
# end
# optimize_pairs(arobust_kernel, arobust_obj, preopt)

# # const candidates = Vector{Float64}[]
# # @time for i in 1:1000
# #     obj, params, ret = NLopt.optimize(opt, Opts.init_vars!(preopt_tracker))
# #     if getfield(NLopt, ret) < 0
# #         continue
# #     end
# #     total_t = preobj(Val((:τ, 0)), params)
# #     dis = preobj(Val((:dis2, 0)), params)
# #     disδ = preobj(Val((:disδ2, 0)), params)
# #     area = 0.0
# #     areaδ = 0.0
# #     for i in 1:7
# #         areai = abs(preobj(Val((:area, i)), params))
# #         areaδi = abs(preobj(Val((:areaδ, i)), params))
# #         area += areai
# #         areaδ += areaδi / areai
# #     end
# #     if abs(dis) < 1e-6 && abs(disδ) < 1e-5 && area >= 100
# #         push!(candidates, params)
# #         @show obj, dis, disδ, area, areaδ, total_t, params[preobj.param.Ωbase]
# #     end
# # end
# # @show length(candidates)

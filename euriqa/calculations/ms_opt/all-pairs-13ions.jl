#!/usr/bin/julia

include("opt-all-pair.jl")
include("ion13-params.jl")

meta = Dict("amp_ratio"=>amp_ratio, "nseg"=>nseg)
candidates = Candidate[]
for d in ARGS
    load_candidates_dir(d, meta=meta, candidates=candidates)
end
println("Loaded $(length(candidates))")
# @assert _meta === nothing || _meta == meta
# const pre_pool = ThreadObjectPool() do
#     return PreOptimizer{nseg}(2π .* fs, 2π .* (fs .+ 0.3); amp_ratio=amp_ratio,
#                               tmin=250, tmax=400, ωmin=2π * 2.28, ωmax=2π * 2.497)
# end
# candidates = @time opt_all_rounds!(pre_pool, nrep, candidates)
# @show length(candidates)
# save_candidates(prefix, candidates, meta)

const mode_info = GateModeInfo(2π .* fs, ηs, bij)
for ion1 in 2:11
    for ion2 in ion1 + 1:12
        weights = get_weights!(mode_info, ion1, ion2)
        @show (ion1, ion2)
        checker = PairChecker(copy(weights), π / 2 / (0.4)^2)
        @show check(checker, candidates)
    end
end

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

#!/usr/bin/julia

mix(s, v) = s ⊻ v
prune(s) = s % 16777216

function next_secret(s)
    s = prune(mix(s, s * 64))
    s = prune(mix(s, s ÷ 32))
    s = prune(mix(s, s * 2048))
    return s
end

function gen_all_banana(s, n)
    ss = [s]
    for i in 1:n
        s = next_secret(s)
        push!(ss, s % 10)
    end
    return ss
end

function find_banana(tgt, ss)
    ns = length(ss)
    nt = length(tgt)
    @inbounds for i in 1:(ns - nt)
        match = true
        for j in 1:nt
            if ss[i + j] - ss[i + j - 1] != tgt[j]
                match = false
                break
            end
        end
        if match
            return ss[i + nt]
        end
    end
    return 0
end

function find_all_tgt(all_ss)
    tgts = Set{NTuple{4,Int}}()
    prices = Dict{NTuple{4,Int},Int}[]
    @inbounds for ss in all_ss
        price = Dict{NTuple{4,Int},Int}()
        for i in 1:(length(ss) - 4)
            tgt = (ss[i + 1] - ss[i], ss[i + 2] - ss[i + 1],
                   ss[i + 3] - ss[i + 2], ss[i + 4] - ss[i + 3])
            push!(tgts, tgt)
            get!(price, tgt, ss[i + 4])
        end
        push!(prices, price)
    end
    @show length(tgts)
    return tgts, prices
end

function sum_banana_for_tgt(tgt, prices)
    return sum(get(price, tgt, 0) for price in prices)
end

function sum_banana(file)
    all_ss = [gen_all_banana(parse(Int, line), 2000) for line in eachline(file)]
    tgts, prices = find_all_tgt(all_ss)
    return maximum(sum_banana_for_tgt(tgt, prices) for tgt in tgts)
end

@show sum_banana(ARGS[1])

# function find_all_tgt(all_ss)
#     tgts = Set{NTuple{4,Int}}()
#     @inbounds for ss in all_ss
#         for i in 1:(length(ss) - 4)
#             push!(tgts, (ss[i + 1] - ss[i],
#                          ss[i + 2] - ss[i + 1],
#                          ss[i + 3] - ss[i + 2],
#                          ss[i + 4] - ss[i + 3]))
#         end
#     end
#     @show length(tgts)
#     return tgts
# end

# function sum_banana_for_tgt(tgt, all_ss)
#     return sum(find_banana(tgt, ss) for ss in all_ss)
# end

# function sum_banana(file)
#     all_ss = [gen_all_banana(parse(Int, line), 2000) for line in eachline(file)]
#     return maximum(sum_banana_for_tgt(tgt, all_ss)
#                    for tgt in find_all_tgt(all_ss))
# end

# @show sum_banana(ARGS[1])

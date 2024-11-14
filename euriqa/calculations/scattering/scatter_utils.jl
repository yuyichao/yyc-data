#!/usr/bin/julia

function rand_sample(rates)
    total_rate = sum(rates)
    t = -log(rand()) / total_rate
    r = total_rate * rand()
    nrates = length(rates)
    for i in 1:nrates
        if r < 0
            return t, i - 1
        end
        r -= rates[i]
    end
    return t, nrates
end

struct ScatterRates
    r0::Float64
    r1::Float64
    rs0::Float64
    rs1::Float64
end

function (rates::ScatterRates)(s0, tmax; maxcount=typemax(Int))
    t = 0.0
    s = s0
    cnt = 0
    while true
        dt, branch = rand_sample(s == 0 ? (rates.r0, rates.rs0) : (rates.r1, rates.rs1))
        if t + dt >= tmax
            return s, cnt
        end
        t += dt
        if branch == 1
            s = 1 - s
        else
            cnt += 1
            if cnt >= maxcount
                return s, cnt
            end
        end
    end
end

mutable struct CountState
    n0::Int
    n1::Int
    CountState() = new(0, 0)
end
function (counter::CountState)((s, cnt))
    if s == 0
        counter.n0 += 1
    else
        counter.n1 += 1
    end
    return
end

struct CountScatter
    counts::Vector{Int}
    CountScatter() = new(Int[])
end
function Base.empty!(counter::CountScatter)
    empty!(counter.counts)
    return counter
end
function (counter::CountScatter)((s, cnt))
    idx = cnt + 1
    cntlen = length(counter.counts)
    if idx > cntlen
        resize!(counter.counts, idx)
        counter.counts[cntlen + 1:idx - 1] .= 0
        counter.counts[idx] = 1
    else
        counter.counts[idx] += 1
    end
    return
end

function repeat_sample(accumulator, n, sampler, args...; kwargs...)
    for i in 1:n
        accumulator(sampler(args...; kwargs...))
    end
end

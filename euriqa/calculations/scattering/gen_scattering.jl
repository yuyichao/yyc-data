#!/usr/bin/julia

using Distributions
using DataStructures
using PyPlot

function total_count(λs)
    sum(rand(Poisson(λ)) for λ in λs)
end

function plot_counts(accum; kws...)
    xs = Int[]
    hs = Int[]
    for (x, h) in accum
        push!(xs, x)
        push!(hs, h)
    end
    bar(xs, hs, width=1; kws...)
end

function uniform_rate(λ0, n)
    return total_count(λ0 / n for i in 1:n)
end
function linear_rate(λ0, n)
    @assert n > 1
    return total_count(2λ0 / n * (i - 1) / (n - 1) for i in 1:n)
end
function binary_rate(λ1, λ2, n)
    if rand(Bool)
        return uniform_rate(λ1, n)
    else
        return uniform_rate(λ2, n)
    end
    # return total_count((i <= n ÷ 2) ? (λ1 / n) : (λ2 / n) for i in 1:n)
end

function random_samples(cb, n)
    accum = Accumulator{Int,Int}()
    for i in 1:n
        push!(accum, cb())
    end
    return accum
end

const single = random_samples(10000) do
    uniform_rate(100, 1)
end
const split = random_samples(10000) do
    uniform_rate(100, 100000)
end
const linear = random_samples(10000) do
    linear_rate(100, 100000)
end
const binary = random_samples(10000) do
    binary_rate(90, 110, 100000)
end

figure()
plot_counts(single, color="C0", alpha=0.7)
plot_counts(split, color="C1", alpha=0.7)
grid()

figure()
plot_counts(linear, color="C0", alpha=0.7)
plot_counts(split, color="C1", alpha=0.7)
plot_counts(binary, color="C2", alpha=0.7)
grid()

show()

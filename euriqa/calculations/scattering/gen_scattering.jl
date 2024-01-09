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

figure()
plot_counts(single, color="C0")
plot_counts(split, color="C1")
grid()

figure()
plot_counts(linear, color="C0")
plot_counts(split, color="C1")
grid()

show()

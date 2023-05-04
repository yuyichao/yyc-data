#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using Statistics
using NaCsData.Fitting: fit_data
using NaCsPlot
using PyPlot

struct Trace
    bin_size::Float64 # in ns
    counts::Vector{Int}
end

read_data(fname) = open(fname) do io
    local bin_sizes
    while !eof(io)
        line = readline(io)
        if line == "#ns/channel"
            line = readline(io)
            bin_sizes = parse.(Float64, split(line))
            break
        end
    end
    @assert @isdefined(bin_sizes)
    found_data = false
    while !eof(io)
        line = readline(io)
        if line == "#counts"
            found_data = true
            break
        end
    end
    @assert found_data
    traces = [Trace(bin_size, Int[]) for bin_size in bin_sizes]
    while !eof(io)
        line = split(readline(io))
        @assert length(line) == length(traces)
        for (trace, word) in zip(traces, line)
            push!(trace.counts, parse(Int, word))
        end
    end
    return traces
end

function truncate_trace!(trace::Trace, freq_MHz)
    nbins = ceil(Int, 1e3 / (trace.bin_size * freq_MHz))
    resize!(trace.counts, nbins)
    return
end

# struct NormalizedTrace
#     bin_size::Float64 # in ns
#     counts::Vector{Float64}
#     uncs::Vector{Float64}
# end

# # For some other measurement we might care about the average power
# # but not for this one...
# function normalize_traces!(traces::AbstractVector{Trace}, bg_trace::Trace)
#     sum_bg = sum(bg_trace.counts)
#     sum_bg_s = sqrt(sum_bg)
#     bg = sum_bg / length(bg_trace.counts)
#     bg_s = sum_bg_s / length(bg_trace.counts)

#     res = NormalizedTrace[]
#     for trace in traces
#         norm_trace = NormalizedTrace(trace.bin_size, Float64[], Float64[])
#         @assert trace.bin_size == bg_trace.bin_size
#         for c in trace.counts
#             push!(norm_trace.counts, c - bg)
#             push!(norm_trace.uncs, sqrt(c + bg_s^2))
#         end
#         m = mean(norm_trace.counts)
#         norm_trace.counts ./= m
#         norm_trace.uncs ./= m
#         push!(res, norm_trace)
#     end
#     return res
# end

# function gen_trace_model(ω)
#     return function (x, p)
#         return 1 .+ p[1] .* sin.(ω .* x) .+ p[2] .* cos.(ω .* x)
#     end
# end

# function fit_trace(trace, freq_MHz)
#     t_ns = (0:length(trace.counts) - 1) .* trace.bin_size
#     model = gen_trace_model(2π * freq_MHz / 1000)
#     fit = fit_data(model, t_ns[2:end - 1], trace.counts[2:end - 1],
#                    trace.uncs[2:end - 1], [0.1, 0])
#     return fit
# end

const data = read_data(joinpath(@__DIR__, "data/20230502.dat"))

const freq_MHz = 118

for trace in data
    truncate_trace!(trace, freq_MHz)
end

const prefix = joinpath(@__DIR__, "imgs", "data_20230502")

function plot_trace(dys, norm_data, fits)
    for (i, dy, trace, fit) in zip(1:length(norm_data), dys, norm_data, fits)
        errorbar((0:length(trace.counts) - 1) .* trace.bin_size, trace.counts,
                 trace.uncs, ls="", color="C$i")
        plot(fit.plotx, fit.ploty, color="C$i", label="$dy")
    end
    grid()
    xlabel("Time (ns)")
    ylabel("Normalized counts")
    legend(fontsize=8, ncol=4)
end

figure()
for (i, trace) in enumerate(data)
    scale = i <= 2 ? 1.08 : 1
    plot((0:length(trace.counts) - 1) .* trace.bin_size, trace.counts .* scale,
         "-", color="C$i")
    # plot(fit.plotx, fit.ploty, color="C$i", label="$dy")
end
grid()
xlabel("Time (ns)")
# ylabel("Normalized counts")
# legend(fontsize=8, ncol=4)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

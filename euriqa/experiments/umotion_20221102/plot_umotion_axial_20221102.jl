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

struct NormalizedTrace
    bin_size::Float64 # in ns
    counts::Vector{Float64}
    uncs::Vector{Float64}
end

# For some other measurement we might care about the average power
# but not for this one...
function normalize_traces!(traces::AbstractVector{Trace}, bg_trace::Trace)
    sum_bg = sum(bg_trace.counts)
    sum_bg_s = sqrt(sum_bg)
    bg = sum_bg / length(bg_trace.counts)
    bg_s = sum_bg_s / length(bg_trace.counts)

    res = NormalizedTrace[]
    for trace in traces
        norm_trace = NormalizedTrace(trace.bin_size, Float64[], Float64[])
        @assert trace.bin_size == bg_trace.bin_size
        for c in trace.counts
            push!(norm_trace.counts, c - bg)
            push!(norm_trace.uncs, sqrt(c + bg_s^2))
        end
        m = mean(norm_trace.counts)
        norm_trace.counts ./= m
        norm_trace.uncs ./= m
        push!(res, norm_trace)
    end
    return res
end

function gen_trace_model(ω)
    return function (x, p)
        return 1 .+ p[1] .* sin.(ω .* x) .+ p[2] .* cos.(ω .* x)
    end
end

function fit_trace(trace, freq_MHz)
    t_ns = (0:length(trace.counts) - 1) .* trace.bin_size
    model = gen_trace_model(2π * freq_MHz / 1000)
    fit = fit_data(model, t_ns[2:end - 1], trace.counts[2:end - 1],
                   trace.uncs[2:end - 1], [0.1, 0])
    return fit
end

const data_1_1 = read_data(joinpath(@__DIR__, "data/Y_-100_20_20_Z_0.dat"))
const data_1_2 = read_data(joinpath(@__DIR__, "data/Y_40_20_100_Z_0.dat"))

const data_2_1 = read_data(joinpath(@__DIR__, "data/Y_-100_20_20_Z_0_2.dat"))
const data_2_2 = read_data(joinpath(@__DIR__, "data/Y_40_20_100_Z_0_2.dat"))

const data_3_1 = read_data(joinpath(@__DIR__, "data/Y_-100_20_20_Z_0_3.dat"))
const data_3_2 = read_data(joinpath(@__DIR__, "data/Y_40_20_100_Z_0_3.dat"))

const data_4_1 = read_data(joinpath(@__DIR__, "data/Y_-100_20_20_Z_0_4.dat"))
const data_4_2 = read_data(joinpath(@__DIR__, "data/Y_40_20_100_Z_0_4.dat"))

const freq_MHz = 46.25
const dys = -100:20:100

for trace in data_1_1
    truncate_trace!(trace, freq_MHz)
end
for trace in data_1_2
    truncate_trace!(trace, freq_MHz)
end

for trace in data_2_1
    truncate_trace!(trace, freq_MHz)
end
for trace in data_2_2
    truncate_trace!(trace, freq_MHz)
end

for trace in data_3_1
    truncate_trace!(trace, freq_MHz)
end
for trace in data_3_2
    truncate_trace!(trace, freq_MHz)
end

for trace in data_4_1
    truncate_trace!(trace, freq_MHz)
end
for trace in data_4_2
    truncate_trace!(trace, freq_MHz)
end

const norm_data_1 = [normalize_traces!(data_1_1[2:end], data_1_1[1]);
                     normalize_traces!(data_1_2[2:end], data_1_2[1])]

const norm_data_2 = [normalize_traces!(data_2_1[2:end], data_2_1[1]);
                     normalize_traces!(data_2_2[2:end], data_2_2[1])]

const norm_data_3 = [normalize_traces!(data_3_1[2:end], data_3_1[1]);
                     normalize_traces!(data_3_2[2:end], data_3_2[1])]

const norm_data_4 = [normalize_traces!(data_4_1[2:end], data_4_1[1]);
                     normalize_traces!(data_4_2[2:end], data_4_2[1])]

const fits_1 = [fit_trace(trace, freq_MHz) for trace in norm_data_1]
const fits_2 = [fit_trace(trace, freq_MHz) for trace in norm_data_2]
const fits_3 = [fit_trace(trace, freq_MHz) for trace in norm_data_3]
const fits_4 = [fit_trace(trace, freq_MHz) for trace in norm_data_4]

const prefix = joinpath(@__DIR__, "imgs", "data_axial_20221102")

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
plot_trace(dys, norm_data_1, fits_1)
NaCsPlot.maybe_save("$(prefix)_0_quantum")

figure()
plot_trace(dys, norm_data_2, fits_2)
NaCsPlot.maybe_save("$(prefix)_0_loading")

figure()
plot_trace(dys, norm_data_3, fits_3)
NaCsPlot.maybe_save("$(prefix)_105_quantum")

figure()
plot_trace(dys, norm_data_4, fits_4)
NaCsPlot.maybe_save("$(prefix)_105_loading")

figure()
errorbar([fit.param[1] for fit in fits_1], [fit.param[2] for fit in fits_1],
         xerr=[fit.unc[1] for fit in fits_1], yerr=[fit.unc[2] for fit in fits_1],
         fmt="o-", label="Quantum")
errorbar([fit.param[1] for fit in fits_2], [fit.param[2] for fit in fits_2],
         xerr=[fit.unc[1] for fit in fits_2], yerr=[fit.unc[2] for fit in fits_2],
         fmt="o-", label="Loading")
grid()
xlim([-0.45, 0.45])
ylim([-0.45, 0.45])
xlabel("Amplitude 1")
ylabel("Amplitude 2")
legend(fontsize=15)
gca()[:set_aspect](1)
NaCsPlot.maybe_save("$(prefix)_amp_phase_0")

figure()
errorbar([fit.param[1] for fit in fits_3], [fit.param[2] for fit in fits_3],
         xerr=[fit.unc[1] for fit in fits_3], yerr=[fit.unc[2] for fit in fits_3],
         fmt="o-", label="Quantum")
errorbar([fit.param[1] for fit in fits_4], [fit.param[2] for fit in fits_4],
         xerr=[fit.unc[1] for fit in fits_4], yerr=[fit.unc[2] for fit in fits_4],
         fmt="o-", label="Loading")
grid()
xlim([-0.45, 0.45])
ylim([-0.45, 0.45])
xlabel("Amplitude 1")
ylabel("Amplitude 2")
legend(fontsize=15)
gca()[:set_aspect](1)
NaCsPlot.maybe_save("$(prefix)_amp_phase_105")

function get_amp_unc(fit)
    amp2 = fit.param[1]^2 + fit.param[2]^2
    # Should use covariance but anyway...
    amp2_s = 2 * abs(fit.param[1]) * fit.unc[1] + 2 * abs(fit.param[1]) * fit.unc[1]
    amp = sqrt(amp2)
    amp_s = amp2_s / amp / 2
    return amp, amp_s
end

function fit_amp(dys, amps)
    minidx = argmin([a[1] for a in amps])
    fit_idxs = (minidx - 2):(minidx + 2)
    function model(x, p)
        return sqrt.(p[1].^2 .+ (p[2] .* (x .- p[3])).^2)
    end
    fit_amps = amps[fit_idxs]
    fit = fit_data(model, dys[fit_idxs], [a[1] for a in fit_amps],
                   [a[2] for a in fit_amps], [0.01, 0.1, dys[minidx]])
    return fit
end

const amp_1 = get_amp_unc.(fits_1)
const amp_2 = get_amp_unc.(fits_2)
const amp_3 = get_amp_unc.(fits_3)
const amp_4 = get_amp_unc.(fits_4)

const fit_amp_1 = fit_amp(dys, amp_1)
const fit_amp_2 = fit_amp(dys, amp_2)
const fit_amp_3 = fit_amp(dys, amp_3)
const fit_amp_4 = fit_amp(dys, amp_4)

figure()
errorbar(dys, [a[1] for a in amp_1], [a[2] for a in amp_1], fmt="C0o", label="Quantum")
plot(fit_amp_1.plotx, fit_amp_1.ploty, "C0-")
errorbar(dys, [a[1] for a in amp_2], [a[2] for a in amp_2], fmt="C1o", label="Loading")
plot(fit_amp_2.plotx, fit_amp_2.ploty, "C1-")
text(-40, 0.4, "\$DV_{min1}=$(fit_amp_1.uncs[3])\$ V/m", color="C0")
text(-40, 0.35, "\$DV_{min2}=$(fit_amp_2.uncs[3])\$ V/m", color="C1")
grid()
legend(fontsize=15, loc="lower left")
xlabel("DY (V/m)")
ylabel("Amplitude")
NaCsPlot.maybe_save("$(prefix)_mindiff_0")

figure()
errorbar(dys, [a[1] for a in amp_3], [a[2] for a in amp_3], fmt="C0o", label="Quantum")
plot(fit_amp_3.plotx, fit_amp_3.ploty, "C0-")
errorbar(dys, [a[1] for a in amp_4], [a[2] for a in amp_4], fmt="C1o", label="Loading")
plot(fit_amp_4.plotx, fit_amp_4.ploty, "C1-")
text(-40, 0.45, "\$DV_{min1}=$(fit_amp_3.uncs[3])\$ V/m", color="C0")
text(-40, 0.4, "\$DV_{min2}=$(fit_amp_4.uncs[3])\$ V/m", color="C1")
grid()
legend(fontsize=15, loc="lower left")
xlabel("DY (V/m)")
ylabel("Amplitude")
NaCsPlot.maybe_save("$(prefix)_mindiff_105")

NaCsPlot.maybe_show()

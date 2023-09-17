#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData.Fitting: fit_data
using Statistics

include("read_data.jl")

function count_cycles(bin)
    data = bin.data
    minv, maxv = extrema(data)
    top_thresh = 0.8 * maxv + 0.2 * minv
    bottom_thresh = 0.2 * maxv + 0.8 * minv

    # 1: raise through top_thresh
    # 2: lower through bottom_thresh
    state = data[1] >= top_thresh ? 2 : 1
    raise_count = 0
    for v in data
        if state == 1
            if v > top_thresh
                raise_count += 1
                state = 2
            end
        else
            if v < bottom_thresh
                state = 1
            end
        end
    end
    spacing = length(data) / raise_count
    period = spacing * bin.dx
    return raise_count, period, Int(maxv) - Int(minv)
end

function find_center_crossings(data)
    minv, maxv = extrema(data)
    mid_v = (minv + maxv) / 2
    top_thresh = 0.8 * maxv + 0.2 * minv
    bottom_thresh = 0.2 * maxv + 0.8 * minv

    crossings = Float64[]
    # 1: raise through top_thresh
    # 2: lower through mid_v
    # 3: lower through bottom_thresh
    # 4: raise through mid_v
    prev_v = data[1]
    state = prev_v >= top_thresh ? 3 : 1
    @inbounds for i in 2:length(data)
        v = Int(data[i])
        if state == 1
            if v > top_thresh
                state = 2
            end
        elseif state == 2
            if v < mid_v
                state = v <= bottom_thresh ? 4 : 3
            end
        elseif state == 3
            if v < bottom_thresh
                state = 4
            end
        else # 4
            if v > mid_v
                state = v >= top_thresh ? 2 : 1
                push!(crossings, i - (v - mid_v) / (v - prev_v))
            end
        end
        prev_v = v
    end
    return crossings
end

# Fit to y = a + b * x
function fit_line(xs, ys)
    m_xy = mean(ys .* xs)
    m_x2 = mean(xs.^2)
    m_x = mean(xs)
    m_y = mean(ys)

    b = (m_xy - m_x * m_y) / (m_x2 - m_x^2)
    a = m_y - b * m_x

    return a, b
end

function estimate_sin2pi_parameter(crossings)
    a, b = fit_line(0:length(crossings) - 1, crossings)
    # a = -ϕ / ω
    # b = 2π / ω
    f = 1 / b
    ϕf = mod(-a * f + 0.5, 1) - 0.5
    return f, ϕf
end

function kernel(t, v0, amp, f, ϕf)
    return v0 .+ amp .* sinpi.(2 .* (f .* t .+ ϕf))
end

function fit_freqs(data, block_size)
    minv, maxv = extrema(data)
    v0_init = (minv + maxv) / 2
    amp_init = (maxv - minv) / 2
    crossings = find_center_crossings(data)
    ncrossings = length(crossings)
    idx_center = Float64[]
    freq = Float64[]
    freq_unc = Float64[]
    function model(x, p)
        return kernel.(x, p[1], p[2], p[3], p[4])
    end

    # last_log = -1

    i = 2
    while true
        c = crossings[i]
        vi1 = floor(Int, c - block_size / 2)
        vi2 = ceil(Int, c + block_size / 2)
        if vi1 < 1
            i += 1
            continue
        elseif vi2 > length(data)
            break
        end
        i1 = i - 1
        while i1 > 1 && c - crossings[i1] <= block_size * 0.55
            i1 -= 1
        end
        i2 = i + 1
        while i2 < ncrossings && crossings[i2] - c <= block_size * 0.55
            i2 += 1
        end
        f_init, ϕf_init = estimate_sin2pi_parameter((@view crossings[i1:i2]) .- c)

        p_init = [v0_init, amp_init, f_init, ϕf_init]
        x = (vi1:vi2) .- c
        v_fit = @view(data[vi1:vi2])
        # if i ÷ 1000 != last_log
        #     last_log = i ÷ 1000
        #     @show i1, i2, vi1, vi2
        # end
        fit = fit_data(model, x, v_fit, p_init, plotx=false)
        push!(idx_center, c)
        push!(freq, fit.param[3])
        push!(freq_unc, fit.unc[3])

        i = i2
    end

    return idx_center, freq, freq_unc
end

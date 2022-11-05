#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using HDF5
using Statistics
using NaCsData.Fitting: fit_data

load_data(fname) = h5open(fname) do io
    return (time=read(io, "time"), value=read(io, "value"),
            tscale=read(io, "tscale"), vscale=read(io, "vscale"))
end

function shrink_block(block, r_start, r_end)
    len = length(block)
    i = round(Int, 1 * (1 - r_start) + len * r_start)
    l = round(Int, 1 * (1 - r_end) + len * r_end)
    return block[i:l]
end

function find_blocks(data)
    threshold = 20
    min_size = 50
    end_threshold = 100

    blocks = UnitRange{Int64}[]
    in_block = false
    block_start = 0
    last_above = 0
    for (i, v) in enumerate(data.value)
        if v < threshold
            if !in_block
                continue
            end
            if i > last_above + end_threshold
                if last_above - block_start > min_size
                    push!(blocks, shrink_block(block_start:last_above, 0.1, 0.9))
                end
                in_block = false
            end
            continue
        end
        last_above = i
        if !in_block
            in_block = true
            block_start = i
        end
    end
    if in_block && last_above - block_start > min_size
        push!(blocks, shrink_block(block_start:last_above, 0.1, 0.9))
    end

    return blocks
end

function find_zeros(ts, data, trigger=0)
    crossings = Float64[]
    for i in 1:length(data) - 1
        d1 = data[i]
        d2 = data[i + 1]
        if d1 <= 0 && d2 > 0
            t1 = ts[i]
            t2 = ts[i + 1]
            push!(crossings, (t1 * d2 - t2 * d1) / (d2 - d1))
        end
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

function estimate_sin2pi_parameter(ts, vs)
    maxv = maximum(vs)
    minv = minimum(vs)
    p1 = (maxv + minv) / 2
    p2 = (maxv - minv) / 2
    crossings = find_zeros(ts, vs, p1)
    a, b = fit_line(0:length(crossings) - 1, crossings)
    # a = -ϕ / ω
    # b = 2π / ω
    ω = 2π / b
    ϕ = mod(-a * ω + π, 2π) - π

    return [p1, p2, ω / 2π, ϕ / 2π]
end

function get_tvblocks(data, blocks)
    tblocks = [data.time[block] for block in blocks]
    vblocks = [data.value[block] for block in blocks]
    return tblocks, vblocks
end

function kernel(t, v0, amp, f, ϕf)
    return v0 .+ amp .* sinpi.(2 .* (f .* t .+ ϕf))
end

function fit_multi_phases(tblocks, vblocks)
    p_ests = estimate_sin2pi_parameter.(tblocks, vblocks)

    p_init = Float64[]
    # first parameter is the global frequency
    push!(p_init, mean([p[1] for p in p_ests]))
    push!(p_init, mean([p[2] for p in p_ests]))
    push!(p_init, mean([p[3] for p in p_ests]))
    for p in p_ests
        push!(p_init, p[4])
    end

    x_real = Tuple{Float64,Int}[]
    for (i, ts) in enumerate(tblocks)
        for t in ts
            push!(x_real, (t, i))
        end
    end

    x_dummy = 1:length(x_real)
    v_fit = reduce(vcat, vec.(vblocks))

    function model_scalar(x, p)
        t, i = x_real[x]
        return kernel(t, p[1], p[2], p[3], p[i + 3])
    end
    model(x, p) = model_scalar.(x, Ref(p))
    fit = fit_data(model, x_dummy, v_fit, p_init, plotx=false)
    return fit
end

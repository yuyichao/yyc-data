#!/usr/bin/julia

function combine_data(reps, vals)
    total_val = zero(eltype(vals))
    total_rep = 0
    for (rep, val) in zip(reps, vals)
        total_val = muladd(rep, val, total_val)
        total_rep += rep
    end
    avg_val = total_val / total_rep
    var = zero(avg_val)
    for (rep, val) in zip(reps, vals)
        var += abs2(rep * (val - avg_val)) / (total_rep - rep)
    end
    return avg_val, sqrt(var / total_rep)
end

struct Collector
    reps::Vector{Int}
    vals::Vector{Float64}
    Collector() = new(Int[], Float64[])
end

function add_data!(c::Collector, rep, v)
    push!(c.reps, rep)
    push!(c.vals, v)
    return
end

function get_data(c::Collector)
    return combine_data(c.reps, c.vals)
end

struct HistCollector
    freq::Dict{Int,Int}
    HistCollector() = new(Dict{Int,Int}())
end

function add_data!(c::HistCollector, _, d)
    for (val, rep) in d
        c.freq[val] = get(c.freq, val, 0) + rep
    end
    return
end

function get_data(c::HistCollector)
    total_val = 0
    total_rep = 0
    for (val, rep) in c.freq
        total_val = muladd(rep, val, total_val)
        total_rep += rep
    end
    avg_val = total_val / total_rep
    var = zero(avg_val)
    for (val, rep) in c.freq
        var += abs2(val - avg_val) * rep
    end
    return avg_val, sqrt(var / (total_rep - 1))
end

struct Line{C}
    collectors::Dict{Int,C}
    xs::Vector{Int}
    vals::Vector{Float64}
    uncs::Vector{Float64}
    Line{C}() where C = new{C}(Dict{Int,C}(), Int[], Float64[], Float64[])
end

function add_point!(l::Line{C}, x, rep, v) where C
    c = get!(C, l.collectors, x)
    add_data!(c, rep, v)
    return
end

function get_data(l::Line)
    if !isempty(l.collectors)
        append!(l.xs, keys(l.collectors))
        sort!(l.xs)
        for s in l.xs
            c = l.collectors[s]
            v, u = get_data(c)
            push!(l.vals, v)
            push!(l.uncs, u)
        end
        empty!(l.collectors)
    end
    return l.xs, l.vals, l.uncs
end

struct LineGroup{C}
    lines::Dict{String,Line{C}}
    LineGroup{C}() where C = new{C}(Dict{String,Line{C}}())
end

function add_test_results!(lg::LineGroup{C}, result, x, rep, conds) where C
    if isa(result, Dict) && first(keys(result)) isa String
        for (k, v) in result
            add_test_results!(lg, v, x, rep, (conds..., k))
        end
    else
        add_point!(get!(Line{C}, lg.lines, join(conds, " ")), x, rep, result)
    end
end

function check_crc32c(result, x, rep, conds)
    if isa(result, Dict) && first(keys(result)) isa String
        for (k, v) in result
            check_crc32c(v, x, rep, (conds..., k))
        end
    else
        if conds == ("ACP_L2", "DDR_WC") || conds == ("ACP_L1", "DDR_WC")
            @assert result != 0
        else
            @assert result == 0
        end
    end
end

function add_results!(all_lines::Dict{<:Any,LineGroup{C}}, results, x, rep, extra_key) where C
    for (k, v) in results
        if k == "crc32c"
            check_crc32c(v, x, rep, ())
            continue
        end
        if extra_key === nothing
            test = k
        else
            test = (k, extra_key)
        end
        lg = get!(LineGroup{C}, all_lines, test)
        add_test_results!(lg, v, x, rep, ())
    end
end

function add_item!(all_lines, item)
    add_results!(all_lines, item["results"], item["size"], item["rep"], get(item, "nbuff", nothing))
end

function add_nbuff_item!(all_lines, item)
    add_results!(all_lines, item["results"], item["nbuff"], item["rep"], nothing)
end

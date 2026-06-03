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

struct Line
    collectors::Dict{Int,Collector}
    sizes::Vector{Int}
    vals::Vector{Float64}
    uncs::Vector{Float64}
    Line() = new(Dict{Int,Collector}(), Int[], Float64[], Float64[])
end

function add_point!(l::Line, size, rep, v)
    c = get!(Collector, l.collectors, size)
    push!(c.reps, rep)
    push!(c.vals, v)
    return
end

function get_data(l::Line)
    if !isempty(l.collectors)
        append!(l.sizes, keys(l.collectors))
        sort!(l.sizes)
        for s in l.sizes
            c = l.collectors[s]
            v, u = combine_data(c.reps, c.vals)
            push!(l.vals, v)
            push!(l.uncs, u)
        end
        empty!(l.collectors)
    end
    return l.sizes, l.vals, l.uncs
end

struct LineGroup
    lines::Dict{String,Line}
    LineGroup() = new(Dict{String,Line}())
end

function add_test_results!(lg, result, size, rep, conds)
    if isa(result, Dict)
        for (k, v) in result
            add_test_results!(lg, v, size, rep, (conds..., k))
        end
    else
        add_point!(get!(Line, lg.lines, join(conds, " ")), size, rep, result)
    end
end

function check_crc32c(result, size, rep, conds)
    if isa(result, Dict)
        for (k, v) in result
            check_crc32c(v, size, rep, (conds..., k))
        end
    else
        if conds == ("ACP_L2", "DDR_WC") || conds == ("ACP_L1", "DDR_WC")
            @assert result != 0
        else
            @assert result == 0
        end
    end
end

function add_results!(all_lines, results, size, rep, nbuff)
    for (k, v) in results
        if k == "crc32c"
            check_crc32c(v, size, rep, ())
            continue
        end
        if nbuff === nothing
            test = k
        else
            test = (k, nbuff)
        end
        lg = get!(LineGroup, all_lines, test)
        add_test_results!(lg, v, size, rep, ())
    end
end

function add_item!(all_lines, item)
    add_results!(all_lines, item["results"], item["size"], item["rep"], get(item, "nbuff", nothing))
end

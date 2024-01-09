#!/usr/bin/julia

mutable struct Disc
    weight::Int
    total_weight::Int
    const holding::Vector{Disc}
    has_imbalance::Bool
    Disc() = new(0, 0, Disc[], false)
end

function scan_disc(disc)
    disc.total_weight = disc.weight
    if isempty(disc.holding)
        return
    end
    for d in disc.holding
        scan_disc(d)
    end
    w = disc.holding[1].total_weight
    for d in disc.holding
        if d.total_weight != w || d.has_imbalance
            disc.has_imbalance = true
        end
        disc.total_weight += d.total_weight
    end
end

function find_wrong_weight(disc, total_weight)
    if !disc.has_imbalance
        return disc.weight + (total_weight - disc.total_weight)
    end
    holding_weight = total_weight - disc.weight
    @assert holding_weight % length(disc.holding) == 0
    w = holding_weight ÷ length(disc.holding)
    res = 0
    found_imbalance = false
    for d in disc.holding
        if w == d.total_weight
            continue
        end
        @assert !found_imbalance
        found_imbalance = true
        res = find_wrong_weight(d, w)
    end
    @assert found_imbalance
    return res
end

function find_wrong_weight(disc)
    if length(disc.holding) == 1
        return find_wrong_weight(disc.holding[1])
    end
    w = 0
    imbalance_count = 0
    for d in disc.holding
        if d.has_imbalance
            imbalance_count += 1
            continue
        end
        if w == 0
            w = d.total_weight
        end
        @assert w == d.total_weight
    end
    @assert imbalance_count == 1
    for d in disc.holding
        if d.has_imbalance
            return find_wrong_weight(d, w)
        end
    end
    error()
end

function balance_tower(file)
    discs = Dict{String,Disc}()
    get_disc(name) = get!(Disc, discs, name)

    subdiscs = Set{Disc}()
    for line in eachline(file)
        line = split(line, " -> ", limit=2)
        m = match(r"([^ ]+) \((\d+)\)", line[1])
        d = get_disc(m[1])
        d.weight = parse(Int, m[2])
        if length(line) > 1
            for name in split(line[2], ",")
                d2 = get_disc(strip(name))
                push!(d.holding, d2)
                push!(subdiscs, d2)
            end
        end
    end

    root = first(setdiff(values(discs), subdiscs))
    scan_disc(root)

    return find_wrong_weight(root)
end

@show balance_tower(ARGS[1])

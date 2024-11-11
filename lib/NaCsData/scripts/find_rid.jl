#!/usr/bin/julia

struct HourData
    min_rid::Int
    max_rid::Int
    data::Dict{Int,String}
end

function load_hour(hourdir)
    data = Dict{Int,String}()
    min_rid = typemax(Int)
    max_rid = typemin(Int)
    for f in readdir(hourdir)
        m = match(r"^([0-9]+)-.*\.h5$", f)
        if m === nothing
            continue
        end
        rid = parse(Int, m[1])
        f = joinpath(hourdir, f)
        isfile(f) || continue
        min_rid = min(min_rid, rid)
        max_rid = max(max_rid, rid)
        data[rid] = f
    end
    if isempty(data)
        return nothing
    end
    return HourData(min_rid, max_rid, data)
end

function find_data(hd::HourData, rid)
    if rid < hd.min_rid
        return -1
    elseif rid > hd.max_rid
        return 1
    end
    file = get(hd.data, rid, nothing)
    if file === nothing
        return 0
    end
    return file
end

struct DateData
    datedir::String
    hours::Vector{String}
    hour_data::Vector{Union{HourData,Nothing}}
end

function load_date(datedir)
    hours = String[]
    for d in readdir(datedir)
        if match(r"^[0-9]*$", d) === nothing
            continue
        end
        isdir(joinpath(datedir, d)) || continue
        push!(hours, d)
    end
    return DateData(datedir, hours,
                    Vector{Union{HourData,Nothing}}(undef, length(hours)))
end

function compute_find_res(lb, ub, nhs, had_valid)
    if !had_valid
        return
    end
    if lb > nhs
        return 1
    elseif ub < 1
        return -1
    end
    return 0
end

function find_data(dd::DateData, rid)
    nhs = length(dd.hours)
    lb = 1
    ub = nhs
    had_valid = false
    while true
        mid0 = (lb + ub) ÷ 2
        mid = mid0
        while true
            if !isassigned(dd.hour_data, mid)
                dd.hour_data[mid] = load_hour(joinpath(dd.datedir, dd.hours[mid]))
            end
            hd = dd.hour_data[mid]
            if hd !== nothing
                res = find_data(hd, rid)
                if isa(res, String)
                    return res
                elseif res == 0
                    return 0
                elseif res < 0
                    ub = mid - 1
                else
                    lb = mid + 1
                end
                had_valid = true
                if lb > ub
                    return compute_find_res(lb, ub, nhs, had_valid)
                end
                break
            end
            if mid >= mid0
                mid += 1
                if mid > ub
                    mid = mid0 - 1
                    ub = mid
                    if mid < lb
                        return compute_find_res(lb, ub, nhs, had_valid)
                    end
                end
            else
                mid -= 1
                if mid < lb
                    lb = ub + 1
                    return compute_find_res(lb, ub, nhs, had_valid)
                end
            end
        end
    end
end

struct AllData
    basedir::String
    dates::Vector{String}
    date_data::Vector{DateData}
end

function AllData(basedir)
    dates = String[]
    for d in readdir(basedir)
        if match(r"^[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]$", d) === nothing
            continue
        end
        isdir(joinpath(basedir, d)) || continue
        push!(dates, d)
    end
    return AllData(basedir, dates, Vector{DateData}(undef, length(dates)))
end

function find_data(ad::AllData, rid)
    nds = length(ad.dates)
    lb = 1
    ub = nds
    while true
        mid0 = (lb + ub) ÷ 2
        mid = mid0
        while true
            if !isassigned(ad.date_data, mid)
                ad.date_data[mid] = load_date(joinpath(ad.basedir, ad.dates[mid]))
            end
            dd = ad.date_data[mid]
            res = find_data(dd, rid)
            if isa(res, String)
                return res
            elseif isa(res, Int)
                if res == 0
                    return
                elseif res < 0
                    ub = mid - 1
                else
                    lb = mid + 1
                end
                if lb > ub
                    return
                end
                break
            end
            if mid >= mid0
                mid += 1
                if mid > ub
                    mid = mid0 - 1
                    ub = mid
                    if mid < lb
                        return
                    end
                end
            else
                mid -= 1
                if mid < lb
                    return
                end
            end
        end
    end
end

const all_data = AllData(ARGS[1])

for rid in ARGS[2:end]
    res = find_data(all_data, parse(Int, rid))
    if res === nothing
        println("$(rid): not found")
    else
        println("$(rid): $(res)")
    end
end

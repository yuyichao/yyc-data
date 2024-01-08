#!/usr/bin/julia

using DataStructures
using Dates

@enum RecordType ShiftStart Sleep Awake

struct Record
    day::Int
    hour::Int
    minute::Int
    rtype::RecordType
    guard_id::Int
end

function load_records(file)
    records = Record[]
    for line in eachline(file)
        m = match(r"\[1518-(\d\d)-(\d\d) (\d\d):(\d\d)\] (wakes up|falls asleep|Guard #(\d+) begins shift)", line)
        month = parse(Int, m[1])
        day = parse(Int, m[2])
        day = Dates.value(Date(1518, month, day) - Date(1518, 1, 1))
        hour = parse(Int, m[3])
        minute = parse(Int, m[4])
        if m[5] == "wakes up"
            record = Record(day, hour, minute, Awake, 0)
        elseif m[5] == "falls asleep"
            record = Record(day, hour, minute, Sleep, 0)
        else
            record = Record(day, hour, minute, ShiftStart, parse(Int, m[6]))
        end
        push!(records, record)
    end
    sort!(records, by=r->(r.day, r.hour, r.minute))
    return records
end

get_night(time) = time[1] + (time[2] > 12)
get_night(r::Record) = get_night((r.day, r.hour, r.minute))

function find_guard(file)
    records = load_records(file)

    guard_record = Dict{Int,Accumulator{Int,Int}}()
    function record_asleep(id, sleep, awake)
        if awake <= sleep
            return
        end
        accum = get!(Accumulator{Int,Int}, guard_record, guard_id)
        for i in sleep:awake - 1
            push!(accum, i)
        end
    end

    prev_night = -1
    sleep_minute = -1
    guard_id = -1
    for record in records
        night = get_night(record)
        if night != prev_night
            if guard_id >= 0 && sleep_minute >= 0
                record_asleep(guard_id, sleep_minute, 60)
            end
            @assert record.rtype == ShiftStart
            prev_night = night
            guard_id = record.guard_id
            sleep_minute = -1
            continue
        end
        @assert record.hour == 0
        if record.rtype == ShiftStart
            if guard_id >= 0 && sleep_minute >= 0
                record_asleep(guard_id, sleep_minute, record.minute)
            end
            guard_id = record.guard_id
            sleep_minute = -1
        elseif record.rtype == Sleep
            sleep_minute = record.minute
        elseif record.rtype == Awake
            @assert sleep_minute >= 0
            record_asleep(guard_id, sleep_minute, record.minute)
            sleep_minute = -1
        end
    end

    guard_id = findmax(sum, guard_record)[2]
    return guard_id * findmax(guard_record[guard_id])[2]
end

@show find_guard(ARGS[1])

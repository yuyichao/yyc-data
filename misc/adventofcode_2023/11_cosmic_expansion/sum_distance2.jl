#!/usr/bin/julia

function expand!(counts)
    gap_start = 1
    offset = 0
    for i in 1:length(counts)
        idx, cnt = counts[i]
        offset += (idx - gap_start) * 999_999
        gap_start = idx + 1
        counts[i] = (idx + offset, cnt)
    end
end

function count_distance(counts, total_count)
    past_points = 0
    last_idx = 0
    distance = 0
    for (idx, cnt) in counts
        distance += (idx - last_idx) * past_points * (total_count - past_points)
        past_points += cnt
        last_idx = idx
    end
    @show total_count, past_points
    return distance
end

function load_image(file)
    galaxy_counts_row = NTuple{2,Int}[]
    galaxy_counts_col = Dict{Int,Int}()
    row = 0
    total_count = 0
    for line in eachline(file)
        row += 1
        row_count = 0
        for i in 1:length(line)
            if line[i] == '#'
                row_count += 1
                galaxy_counts_col[i] = get(galaxy_counts_col, i, 0) + 1
            end
        end
        if row_count != 0
            push!(galaxy_counts_row, (row, row_count))
        end
        total_count += row_count
    end

    galaxy_counts_col = sort!([(i, c) for (i, c) in galaxy_counts_col])
    expand!(galaxy_counts_col)
    expand!(galaxy_counts_row)

    return (count_distance(galaxy_counts_col, total_count) +
        count_distance(galaxy_counts_row, total_count))
end

@show load_image(ARGS[1])

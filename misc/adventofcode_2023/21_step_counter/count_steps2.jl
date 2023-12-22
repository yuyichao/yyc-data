#!/usr/bin/julia

function propagate!(garden_map, x0, y0, frontier)
    step = garden_map[y0, x0]
    next_step = step + 1
    for (x, y) in ((x0 + 1, y0), (x0 - 1, y0), (x0, y0 + 1), (x0, y0 - 1))
        if !checkbounds(Bool, garden_map, y, x) || garden_map[y, x] <= next_step
            continue
        end
        garden_map[y, x] = next_step
        push!(frontier, (x, y))
    end
end

function propagate_all!(garden_map, frontier)
    while !isempty(frontier)
        (x, y) = popfirst!(frontier)
        propagate!(garden_map, x, y, frontier)
    end
end

function reset_map!(garden_map)
    for idx in eachindex(garden_map)
        if garden_map[idx] != -1
            garden_map[idx] = typemax(Int)
        end
    end
end

struct PropagateResult
    step_counts::Vector{Int} # Starting from 0
    function PropagateResult(garden_map)
        step_counts = Int[]
        for v in garden_map
            if v == -1 || v == typemax(Int)
                continue
            end
            idx = v + 1
            ncounts = length(step_counts)
            if idx > ncounts
                resize!(step_counts, idx)
                step_counts[ncounts + 1:idx - 1] .= 0
                step_counts[idx] = 1
            else
                step_counts[idx] += 1
            end
        end
        return new(step_counts)
    end
end

function count_result(result::PropagateResult, max_step)
    s = 0
    for (i, c) in enumerate(result.step_counts)
        x = i - 1
        if x <= max_step && (max_step - x) % 2 == 0
            s += c
        end
    end
    return s
end
get_max_step(result::PropagateResult) = length(result.step_counts) - 1

function count_straight(block_result, max_steps, side_len)
    blk_max_step = get_max_step(block_result)
    nfull = (max_steps - blk_max_step) ÷ side_len + 1
    s = 0
    if nfull > 0
        nodd = nfull ÷ 2
        neven = nfull - nodd
        s += (neven * count_result(block_result, max_steps) +
            nodd * count_result(block_result, max_steps - side_len))
        max_steps -= nfull * side_len
    end
    while max_steps >= 0
        s += count_result(block_result, max_steps)
        max_steps -= side_len
    end
    return s
end

function count_corner(block_result, max_steps, ncols, nrows)
    s = 0
    while max_steps >= 0
        s += count_straight(block_result, max_steps, ncols)
        max_steps -= nrows
    end
    return s
end

function count_steps(garden_map::Matrix, start_pos, max_steps)
    for i in 1:size(garden_map, 1)
        @assert garden_map[i, 1] != -1
        @assert garden_map[i, start_pos[1]] != -1
        @assert garden_map[i, end] != -1
    end
    for i in 1:size(garden_map, 2)
        @assert garden_map[1, i] != -1
        @assert garden_map[start_pos[2], i] != -1
        @assert garden_map[end, i] != -1
    end

    nrows, ncols = size(garden_map)

    xs = (1, start_pos[1], ncols)
    ys = (1, start_pos[2], nrows)

    # @assert xs[2] - xs[1] == xs[3] - xs[2]
    # @assert ys[2] - ys[1] == ys[3] - ys[2]

    frontier = NTuple{2,Int}[]

    corner_results = [begin
                          push!(frontier, (y0, x0))
                          reset_map!(garden_map)
                          garden_map[x0, y0] = 0
                          propagate_all!(garden_map, frontier)
                          PropagateResult(garden_map)
                      end for y0 in ys, x0 in xs]

    s = count_result(corner_results[2, 2], max_steps)

    s += count_straight(corner_results[1, 2], max_steps - start_pos[2], nrows)
    s += count_straight(corner_results[2, 1], max_steps - start_pos[1], ncols)
    s += count_straight(corner_results[3, 2],
                        max_steps - (nrows - start_pos[2] + 1), nrows)
    s += count_straight(corner_results[2, 3],
                        max_steps - (ncols - start_pos[1] + 1), ncols)

    s += count_corner(corner_results[1, 1], max_steps - start_pos[2] - start_pos[1],
                      ncols, nrows)
    s += count_corner(corner_results[1, 3],
                      max_steps - start_pos[2] - (ncols - start_pos[1] + 1),
                      ncols, nrows)
    s += count_corner(corner_results[3, 1],
                      max_steps - start_pos[1] - (nrows - start_pos[2] + 1),
                      ncols, nrows)
    s += count_corner(corner_results[3, 3],
                      max_steps - (ncols - start_pos[1] + 1) -
                          (nrows - start_pos[2] + 1),
                      ncols, nrows)

    return s
end

function count_steps(file)
    lines = readlines(file)
    M = Matrix{Int}(undef, length(lines), length(lines[1]))
    start_pos = (0, 0)
    for i in 1:length(lines)
        line = lines[i]
        for j in 1:length(line)
            c = line[j]
            if c == '#'
                M[i, j] = -1
            elseif c == 'S'
                M[i, j] = 0
                start_pos = (j, i)
            else
                @assert c == '.'
                M[i, j] = typemax(Int)
            end
        end
    end

    max_steps = 26501365
    count_steps(M, start_pos, max_steps)

    # corner_results = PropagateResult[]

    # # 12
    # # 34
    # for (y, x) in ((1, 1), (1, size(M, 2)),
    #                (size(M, 1), 1), size(M))
    #     reset_map!(M)
    #     M[y, x] = 0
    #     push!(frontier, (x, y))
    #     propagate_all!(M, frontier)
    #     push!(corner_results, PropagateResult(M))
    # end

    # @show corner_results

    # display(M)
    # return count(x->x <= 64 && x % 2 == 0, M)
end

@show count_steps(ARGS[1])

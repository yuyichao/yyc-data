#!/usr/bin/julia

function calc_path(file)
    lines = readlines(file)
    nrows = length(lines)
    ncols = length(lines[1])
    M = zeros(Int, nrows, ncols)
    start_pos = (0, 0)
    end_pos = (0, 0)
    for row in 1:nrows
        line = lines[row]
        for col in 1:ncols
            c = line[col]
            if c == '.'
            elseif c == '#'
                M[row, col] = 1
            elseif c == 'S'
                start_pos = (row, col)
            elseif c == 'E'
                end_pos = (row, col)
            else
                error("....")
            end
        end
    end
    scores = (fill(-1, nrows, ncols),
              fill(-1, nrows, ncols),
              fill(-1, nrows, ncols),
              fill(-1, nrows, ncols))
    dirs = ((0, 1), (0, -1), (1, 0), (-1, 0))

    todo = Set(((start_pos..., 1),))
    scores[1][start_pos...] = 0
    while !isempty(todo)
        row, col, dir = pop!(todo)
        pos = row, col
        score = scores[dir][row, col]
        for i in 1:4
            if i == dir
                new_pos = pos .+ dirs[dir]
                new_dir = dir
                if M[new_pos...] != 0
                    continue
                end
                new_score = score + 1
            else
                new_pos = pos
                new_dir = i
                new_score = score + 1000
            end
            if (scores[new_dir][new_pos...] == -1 ||
                scores[new_dir][new_pos...] > new_score)
                scores[new_dir][new_pos...] = new_score
                push!(todo, (new_pos..., new_dir))
            end
        end
    end

    min_score = min(scores[1][end_pos...], scores[2][end_pos...],
                    scores[3][end_pos...], scores[4][end_pos...])

    empty!(todo)
    for i in 1:4
        if scores[i][end_pos...] == min_score
            push!(todo, (end_pos..., i))
        end
    end
    path_pos = Set((end_pos,))
    while !isempty(todo)
        row, col, dir = pop!(todo)
        pos = row, col
        score = scores[dir][row, col]
        for i in 1:4
            if i == dir
                prev_pos = pos .- dirs[dir]
                prev_dir = dir
                if M[prev_pos...] != 0
                    continue
                end
                prev_score = score - 1
            else
                prev_pos = pos
                prev_dir = i
                prev_score = score - 1000
            end
            if scores[prev_dir][prev_pos...] == prev_score
                push!(path_pos, prev_pos)
                push!(todo, (prev_pos..., prev_dir))
            end
        end
    end
    @assert start_pos in path_pos
    return length(path_pos)
end

@show calc_path(ARGS[1])

#!/usr/bin/julia

function min_score(M, start_pos)
    nrows, ncols = size(M)
    scores = fill(-1, nrows, ncols)

    todo = Set((start_pos,))
    scores[start_pos...] = 0
    while !isempty(todo)
        pos = pop!(todo)
        score = scores[pos...]
        for dp in ((0, 1), (0, -1), (1, 0), (-1, 0))
            new_pos = pos .+ dp
            new_score = score + 1
            if M[new_pos...] == 0 && (scores[new_pos...] == -1 ||
                scores[new_pos...] > new_score)
                scores[new_pos...] = new_score
                push!(todo, new_pos)
            end
        end
    end
    return scores
end

function count_cheat(file)
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
    forward_scores = min_score(M, start_pos)
    backward_scores = min_score(M, end_pos)
    old_score = forward_scores[end_pos...]
    @assert old_score == backward_scores[start_pos...]

    s = 0
    for row1 in 2:nrows - 1
        for col1 in 2:ncols - 1
            if M[row1, col1] == 1
                continue
            end
            for row2 in 2:nrows - 1
                for col2 in 2:ncols - 1
                    if M[row2, col2] == 1 || (row1 == row2 && col1 == col2)
                        continue
                    end
                    d = abs(row1 - row2) + abs(col1 - col2)
                    if d > 20
                        continue
                    end
                    new_score = (forward_scores[row1, col1] + d +
                        backward_scores[row2, col2])
                    if new_score <= old_score - 100
                        s += 1
                    end
                end
            end
        end
    end
    return s
end

@show count_cheat(ARGS[1])

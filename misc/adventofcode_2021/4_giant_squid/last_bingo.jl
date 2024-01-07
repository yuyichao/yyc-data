#!/usr/bin/julia

function has_bingo(marked)
    for i in 1:5
        all_marked = true
        for j in 1:5
            if !marked[i, j]
                all_marked = false
                break
            end
        end
        if all_marked
            return true
        end
        all_marked = true
        for j in 1:5
            if !marked[j, i]
                all_marked = false
                break
            end
        end
        if all_marked
            return true
        end
    end
    return false
end

function mark_number(board, marked, num)
    for i in 1:5
        for j in 1:5
            if board[j, i] == num
                marked[j, i] = true
            end
        end
    end
end

function find_bingo(board, numbers)
    @assert size(board) == (5, 5)
    marked = zeros(Bool, 5, 5)
    for (i, num) in enumerate(numbers)
        mark_number(board, marked, num)
        if has_bingo(marked)
            return (i, num * sum(v for (v, m) in zip(board, marked) if !m))
        end
    end
    return typemax(Int), 0
end

function first_bingo(file)
    lines = eachline(file)
    numbers = parse.(Int, split(first(lines), ','))
    board = Matrix{Int}(undef, 5, 5)
    winning_round = 0
    winning_score = 0
    for line in lines
        @assert isempty(line)
        board[1, :] = parse.(Int, split(first(lines)))
        board[2, :] = parse.(Int, split(first(lines)))
        board[3, :] = parse.(Int, split(first(lines)))
        board[4, :] = parse.(Int, split(first(lines)))
        board[5, :] = parse.(Int, split(first(lines)))
        r, s = find_bingo(board, numbers)
        if r > winning_round
            winning_round = r
            winning_score = s
        end
    end
    return winning_score
end

@show first_bingo(ARGS[1])

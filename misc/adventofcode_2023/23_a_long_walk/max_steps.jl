#!/usr/bin/julia

const LEFT = 1
const RIGHT = 2
const UP = 3
const DOWN = 4

function check_pos(M, visited, x, y)
    if !checkbounds(Bool, M, y, x) || M[y, x] == -1
        return false
    end
    if (x, y) in visited
        return false
    end
    return true
end

function trace_max_steps(M, pos, step, end_pos, visited, max_steps)
    x, y = pos

    next_pos = NTuple{2,Int}[]

    while true
        if max_steps[y, x] < step
            max_steps[y, x] = step
        end
        if end_pos == (x, y)
            return
        end
        push!(visited, (x, y))
        v = M[y, x]
        step += 1
        if v == LEFT
            x = x - 1
            if !check_pos(M, visited, x, y)
                return
            end
        elseif v == RIGHT
            x = x + 1
            if !check_pos(M, visited, x, y)
                return
            end
        elseif v == UP
            y = y - 1
            if !check_pos(M, visited, x, y)
                return
            end
        elseif v == DOWN
            y = y + 1
            if !check_pos(M, visited, x, y)
                return
            end
        else
            @assert v == 0
            empty!(next_pos)
            for (x2, y2) in ((x - 1, y), (x + 1, y),
                             (x, y - 1), (x, y + 1))
                if !check_pos(M, visited, x2, y2)
                    continue
                end
                push!(next_pos, (x2, y2))
            end
            if isempty(next_pos)
                return
            end
            if length(next_pos) > 1
                for i in 2:length(next_pos)
                    trace_max_steps(M, next_pos[i], step, end_pos,
                                    copy(visited), max_steps)
                end
            end
            x, y = next_pos[1]
        end
    end
end

function trace_max_steps(M, start_pos, end_pos)
    nrows, ncols = size(M)
    max_steps = zeros(Int, nrows, ncols)
    visited = Set{NTuple{2,Int}}()
    trace_max_steps(M, start_pos, 0, end_pos, visited, max_steps)
    return max_steps[end_pos[2], end_pos[1]]
end

function max_steps(file)
    frontier = NTuple{2,Int}[]

    lines = readlines(file)
    nrows = length(lines)
    ncols = length(lines[1])

    M = Matrix{Int8}(undef, nrows, ncols)
    start_pos = (0, 0)
    end_pos = (0, 0)
    for i in 1:nrows
        line = lines[i]
        for j in 1:ncols
            c = line[j]
            if c == '#'
                M[i, j] = -1
            elseif c == '<'
                M[i, j] = LEFT
            elseif c == '>'
                M[i, j] = RIGHT
            elseif c == '^'
                M[i, j] = UP
            elseif c == 'v'
                M[i, j] = DOWN
            else
                @assert c == '.'
                M[i, j] = 0
                if i == 1
                    @assert start_pos == (0, 0)
                    start_pos = (j, i)
                elseif i == nrows
                    @assert end_pos == (0, 0)
                    end_pos = (j, i)
                end
            end
        end
    end

    return trace_max_steps(M, start_pos, end_pos)
end

@show max_steps(ARGS[1])

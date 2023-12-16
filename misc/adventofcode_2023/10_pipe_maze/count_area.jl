#!/usr/bin/julia

const dir_map = Dict('|'=>Dict(3=>3, 4=>4),
                     '-'=>Dict(1=>1, 2=>2),
                     'L'=>Dict(4=>2, 1=>3),
                     'J'=>Dict(4=>1, 2=>3),
                     '7'=>Dict(3=>1, 2=>4),
                     'F'=>Dict(3=>2, 1=>4))

# direction: left, right, top, down
function check_next(lines, x, y, dir)
    nlines = length(lines)
    if y > nlines || y < 1
        return
    end
    line = lines[y]
    ncol = length(line)
    if x > ncol || x < 1
        return
    end

    site = line[x]
    if site == 'S'
        return true
    elseif site == '.'
        return
    end
    return get(dir_map[site], dir, nothing)
end

function track_path(lines, x, y, dir)
    path = NTuple{2,Int}[]
    while true
        if dir == 1
            x, y = x - 1, y
        elseif dir == 2
            x, y = x + 1, y
        elseif dir == 3
            x, y = x, y - 1
        elseif dir == 4
            x, y = x, y + 1
        else
            error("Invalid direction $dir")
        end
        push!(path, (x, y))
        dir = check_next(lines, x, y, dir)
        if dir === true
            return path
        elseif dir === nothing
            return
        end
    end
end

function find_dir(xy1, xy2)
    x1, y1 = xy1
    x2, y2 = xy2
    if x1 < x2
        return 2
    elseif x2 < x1
        return 1
    elseif y1 < y2
        return 4
    else
        return 3
    end
end

# state:
#   0: outside
#   1: outside, below line
#   2: outside, above line
#   3: inside
function count_line_area(raw_line, line)
    area = 0

    state = 0
    start_x = 0

    sort!(line)
    for (x, typ) in line
        if typ == '-'
            @assert state == 1 || state == 2
        elseif typ == '|'
            @assert state == 0 || state == 3
            if state == 0
                state = 3
                start_x = x + 1
            else
                state = 0
                area += x - start_x
            end
        elseif typ == 'L'
            @assert state == 0 || state == 3
            if state == 0
                state = 1
            else
                state = 2
                area += x - start_x
            end
        elseif typ == 'F'
            @assert state == 0 || state == 3
            if state == 0
                state = 2
            else
                state = 1
                area += x - start_x
            end
        elseif typ == 'J'
            @assert state == 1 || state == 2
            if state == 2
                state = 3
                start_x = x + 1
            else
                state = 0
            end
        elseif typ == '7'
            @assert state == 1 || state == 2
            if state == 1
                state = 3
                start_x = x + 1
            else
                state = 0
            end
        else
            error("Unknown type $typ")
        end
    end
    return area
end

function trace_path(lines, path)
    edges = Dict{Int,Vector{Tuple{Int,Char}}}()
    for i in 1:length(path) - 1
        x, y = path[i]
        push!(get!(Vector{Tuple{Int,Char}}, edges, y), (x, lines[y][x]))
    end
    start_dir = find_dir(path[end], path[1])
    end_dir = find_dir(path[end - 1], path[end])

    if (end_dir, start_dir) == (1, 1) ||  (end_dir, start_dir) == (2, 2)
        start_edge = '-'
    elseif (end_dir, start_dir) == (3, 3) ||  (end_dir, start_dir) == (4, 4)
        start_edge = '|'
    elseif (end_dir, start_dir) == (1, 3) ||  (end_dir, start_dir) == (4, 2)
        start_edge = 'L'
    elseif (end_dir, start_dir) == (1, 4) ||  (end_dir, start_dir) == (3, 2)
        start_edge = 'F'
    elseif (end_dir, start_dir) == (2, 3) ||  (end_dir, start_dir) == (4, 1)
        start_edge = 'J'
    elseif (end_dir, start_dir) == (2, 4) ||  (end_dir, start_dir) == (3, 1)
        start_edge = '7'
    end
    x, y = path[end]
    push!(get!(Vector{Tuple{Int,Char}}, edges, y), (x, start_edge))

    return sum(count_line_area(lines[lineno], line) for (lineno, line) in edges)
end

function count_area(file)
    lines = readlines(file)
    nlines = length(lines)
    sloc = (0, 0)
    for i in 1:nlines
        line = lines[i]
        j = findfirst('S', line)
        if j !== nothing
            sloc = (j, i)
            break
        end
    end

    for i in 1:4
        path = track_path(lines, sloc..., i)
        if path === nothing
            continue
        end
        return trace_path(lines, path)
    end
end

@show count_area(ARGS[1])

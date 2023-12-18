#!/usr/bin/julia

const corner_map = Dict(('R', 'R')=>'-',
                        ('L', 'L')=>'-',
                        ('U', 'U')=>'|',
                        ('D', 'D')=>'|',
                        ('R', 'U')=>'J',
                        ('D', 'L')=>'J',
                        ('R', 'D')=>'7',
                        ('U', 'L')=>'7',
                        ('L', 'U')=>'L',
                        ('D', 'R')=>'L',
                        ('L', 'D')=>'F',
                        ('U', 'R')=>'F')

function trace_edge(instructions)
    last_dir = instructions[end][1]
    x = 0
    y = 0
    edges = Dict{Int,Vector{Tuple{Int,Char}}}()
    for (dir, steps) in instructions
        push!(get!(Vector{Tuple{Int,Char}}, edges, y),
              (x, corner_map[(last_dir, dir)]))
        last_dir = dir
        if dir == 'R'
            dx = 1
            dy = 0
            edge = '-'
        elseif dir == 'L'
            dx = -1
            dy = 0
            edge = '-'
        elseif dir == 'U'
            dx = 0
            dy = -1
            edge = '|'
        elseif dir == 'D'
            dx = 0
            dy = 1
            edge = '|'
        else
            error()
        end
        for i in 1:steps - 1
            push!(get!(Vector{Tuple{Int,Char}}, edges, y + dy * i),
                  (x + dx * i, edge))
        end
        x += dx * steps
        y += dy * steps
    end
    return edges
end

# state:
#   0: outside
#   1: inside, below line
#   2: inside, above line
#   3: inside
function count_line_area(line)
    area = 0

    sort!(line)

    state = 0
    start_x = line[1][1] - 1

    for (x, typ) in line
        if typ == '-'
            @assert state == 1 || state == 2
            area += 1
        elseif typ == '|'
            @assert state == 0 || state == 3
            if state == 0
                state = 3
                start_x = x
            else
                state = 0
                area += x - start_x + 1
            end
        elseif typ == 'L'
            @assert state == 0 || state == 3
            if state == 0
                state = 2
                area += 1
            else
                state = 1
                area += x - start_x + 1
            end
        elseif typ == 'F'
            @assert state == 0 || state == 3
            if state == 0
                state = 1
                area += 1
            else
                state = 2
                area += x - start_x + 1
            end
        elseif typ == 'J'
            @assert state == 1 || state == 2
            if state == 1
                state = 3
                start_x = x
            else
                state = 0
                area += 1
            end
        elseif typ == '7'
            @assert state == 1 || state == 2
            if state == 2
                state = 3
                start_x = x
            else
                state = 0
                area += 1
            end
        else
            error("Unknown type $typ")
        end
    end
    return area
end

function trace_path(instructions)
    edges = trace_edge(instructions)
    return sum(count_line_area(line) for (lineno, line) in edges)
end

function count_area(file)
    instructions = Tuple{Char,Int}[]
    for line in eachline(file)
        m = match(r"^([RLUD]) ([0-9]+) ", line)
        push!(instructions, (m[1][1], parse(Int, m[2])))
    end
    return trace_path(instructions)
end

@show count_area(ARGS[1])

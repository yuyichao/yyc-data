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
    count = 0
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
        count += 1
        dir = check_next(lines, x, y, dir)
        if dir === true
            return count
        elseif dir === nothing
            return 0
        end
    end
end

function count_path(file)
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
    return [track_path(lines, sloc..., i) for i in 1:4]
end

@show count_path(ARGS[1])

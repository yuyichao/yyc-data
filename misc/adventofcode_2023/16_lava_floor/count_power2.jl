#!/usr/bin/julia

const LEFT = 0x1
const RIGHT = 0x2
const UP = 0x4
const DOWN = 0x8

macro check_dir(dir, left, right, up, down)
    quote
        _dir = $(esc(dir))
        if _dir == LEFT
            $(esc(left))
        elseif _dir == RIGHT
            $(esc(right))
        elseif _dir == UP
            $(esc(up))
        elseif _dir == DOWN
            $(esc(down))
        else
            error(_dir)
        end
    end
end

next_pixel(x, y, dir) =
    @check_dir(dir, (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))

function check_and_record(pixel_state, x, y, dir)
    cur_pixel = pixel_state[y, x]
    if (cur_pixel & dir) != 0
        return false
    end
    pixel_state[y, x] = cur_pixel | dir
    return true
end

@inline function handle_new_pixel(lines, pixel_state, x, y, dir)
    if x < 1 || y < 1 || x > size(pixel_state, 2) || y > size(pixel_state, 1)
        return false
    end
    if !check_and_record(pixel_state, x, y, dir[])
        return false
    end
    p = lines[y][x]
    if p == '.'
        return true
    elseif p == '/'
        dir[] = @check_dir(dir[], DOWN, UP, RIGHT, LEFT)
    elseif p == '\\'
        dir[] = @check_dir(dir[], UP, DOWN, LEFT, RIGHT)
    elseif p == '-'
        if dir[] == UP || dir[] == DOWN
            if check_and_record(pixel_state, x, y, LEFT)
                propagate!(lines, pixel_state, x, y, LEFT)
            end
            if !check_and_record(pixel_state, x, y, RIGHT)
                return false
            end
            dir[] = RIGHT
        end
    elseif p == '|'
        if dir[] == LEFT || dir[] == RIGHT
            if check_and_record(pixel_state, x, y, UP)
                propagate!(lines, pixel_state, x, y, UP)
            end
            if !check_and_record(pixel_state, x, y, DOWN)
                return false
            end
            dir[] = DOWN
        end
    else
        error(p)
    end
    return true
end

function propagate!(lines, pixel_state, x, y, dir)
    dir = Ref(dir)
    while true
        x, y = next_pixel(x, y, dir[])
        if !handle_new_pixel(lines, pixel_state, x, y, dir)
            return
        end
    end
end

function count_power(file)
    lines = readlines(file)
    nrows = length(lines)
    ncols = length(lines[1])
    pixel_state = Array{UInt8}(undef, nrows, ncols)

    max_power = 0

    dir = Ref(DOWN)

    @time for x in 1:ncols
        pixel_state .= 0x0
        dir[] = DOWN
        if handle_new_pixel(lines, pixel_state, x, 1, dir)
            propagate!(lines, pixel_state, x, 1, dir[])
        end
        max_power = max(max_power, count(!=(0), pixel_state))

        pixel_state .= 0x0
        dir[] = UP
        if handle_new_pixel(lines, pixel_state, x, nrows, dir)
            propagate!(lines, pixel_state, x, nrows, dir[])
        end
        max_power = max(max_power, count(!=(0), pixel_state))
    end

    @time for y in 1:nrows
        pixel_state .= 0x0
        dir[] = RIGHT
        if handle_new_pixel(lines, pixel_state, 1, y, dir)
            propagate!(lines, pixel_state, 1, y, dir[])
        end
        max_power = max(max_power, count(!=(0), pixel_state))

        pixel_state .= 0x0
        dir[] = LEFT
        if handle_new_pixel(lines, pixel_state, ncols, y, dir)
            propagate!(lines, pixel_state, ncols, y, dir[])
        end
        max_power = max(max_power, count(!=(0), pixel_state))
    end

    return max_power
end

@show count_power(ARGS[1])

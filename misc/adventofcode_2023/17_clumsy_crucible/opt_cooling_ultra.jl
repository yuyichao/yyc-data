#!/usr/bin/julia

const LEFT = 1
const RIGHT = 2
const UP = 3
const DOWN = 4

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
next_site(x, y, dir) =
    @check_dir(dir, (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
next_turn_dir(dir) =
    @check_dir(dir, (UP, DOWN), (UP, DOWN), (LEFT, RIGHT), (LEFT, RIGHT))

struct Site
    cool::Int8
    # The three numbers in each column is the minimum cooling we can achieve
    # when we are still allowed to move in that direction out from this site
    # by 1-10 steps.
    min_cooling::Matrix{Int}
    function Site(cool)
        return new(cool, fill(typemax(Int), 10, 4))
    end
end

function load_lines(lines)
    nrows = length(lines)
    ncols = length(lines[1])
    sites = [Site(parse(Int, lines[i][j])) for i in 1:nrows, j in 1:ncols]
    sites[1, 1].min_cooling[4:10, RIGHT] .= 0
    sites[1, 1].min_cooling[4:10, DOWN] .= 0
    return sites
end

function check_site_dir(sites, workset, x, y, dir)
    site = sites[y, x]
    x2, y2 = next_site(x, y, dir)
    if x2 < 1 || y2 < 1 || x2 > size(sites, 2) || y2 > size(sites, 1)
        return
    end
    site2 = sites[y2, x2]
    changed = false
    for i in 1:9
        cool = site.min_cooling[i + 1, dir]
        if cool == typemax(Int)
            continue
        end
        new_cool = cool + site2.cool
        if new_cool < site2.min_cooling[i, dir]
            changed = true
            site2.min_cooling[i, dir] = new_cool
        end
    end
    cool = site.min_cooling[1, dir]
    if cool < typemax(Int)
        new_cool = cool + site2.cool
        for dir2 in next_turn_dir(dir)
            for i in 4:10
                if new_cool < site2.min_cooling[i, dir2]
                    changed = true
                    site2.min_cooling[i, dir2] = new_cool
                end
            end
        end
    end
    if changed
        push!(workset, (x2, y2))
    end
end

function check_site(sites, workset, x, y)
    check_site_dir(sites, workset, x, y, LEFT)
    check_site_dir(sites, workset, x, y, RIGHT)
    check_site_dir(sites, workset, x, y, UP)
    check_site_dir(sites, workset, x, y, DOWN)
end

function propagate(sites, workset)
    while !isempty(workset)
        (x, y) = pop!(workset)
        check_site(sites, workset, x, y)
    end
end

function opt_cooling(file)
    sites = load_lines(readlines(file))
    propagate(sites, Set(((1, 1),)))
    return minimum(sites[end, end].min_cooling)
end
@show opt_cooling(ARGS[1])

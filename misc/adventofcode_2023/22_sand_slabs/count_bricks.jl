#!/usr/bin/julia

struct Brick
    r0::NTuple{3,Int}
    r1::NTuple{3,Int}
    Brick(r0, r1) = new(min.(r0, r1), max.(r0, r1))
end

function read_bricks(file)
    bricks = Brick[]
    for line in eachline(file)
        m = match(r"(\d+),(\d+),(\d+)~(\d+),(\d+),(\d+)", line)
        push!(bricks, Brick((parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])),
                            (parse(Int, m[4]), parse(Int, m[5]), parse(Int, m[6]))))
    end
    return bricks
end

range_overlap(rg1, rg2) = !(rg1[1] > rg2[2] || rg2[1] > rg1[2])

xy_overlap(brick1, brick2) =
    (range_overlap((brick1.r0[1], brick1.r1[1]), (brick2.r0[1], brick2.r1[1])) &&
    range_overlap((brick1.r0[2], brick1.r1[2]), (brick2.r0[2], brick2.r1[2])))

function drop_bricks!(bricks)
    nbricks = length(bricks)
    for i in 1:nbricks
        brick = bricks[i]
        min_z = 1
        min_new_idx = 1
        for j in i - 1:-1:1
            brick2 = bricks[j]
            if xy_overlap(brick, brick2)
                min_new_idx = j + 1
                min_z = brick2.r1[3] + 1
                break
            end
        end
        @assert min_z <= brick.r0[3]
        if min_z == brick.r0[3]
            continue
        end
        brick = Brick((brick.r0[1], brick.r0[2], min_z),
                      (brick.r1[1], brick.r1[2], min_z - brick.r0[3] + brick.r1[3]))
        if min_new_idx == i
            bricks[i] = brick
            continue
        end
        new_i = (min_new_idx - 1 +
            searchsortedfirst(@view(bricks[min_new_idx:i - 1]), brick,
                              by=b->(b.r1[3], b.r0[3])))
        for j in i - 1:-1:new_i
            bricks[j + 1] = bricks[j]
        end
        bricks[new_i] = brick
    end
end

function compute_breakable(bricks)
    nbricks = length(bricks)
    breakable = ones(Bool, nbricks)
    for i in 1:nbricks
        brick = bricks[i]
        unique_j = 0
        for j in i - 1:-1:1
            brick2 = bricks[j]
            if brick2.r1[3] >= brick.r0[3]
                continue
            elseif brick2.r1[3] < brick.r0[3] - 1
                break
            end
            @assert brick2.r1[3] == brick.r0[3] - 1
            if xy_overlap(brick, brick2)
                if unique_j != 0
                    unique_j = 0
                    break
                end
                unique_j = j
            end
        end
        if unique_j != 0
            breakable[unique_j] = false
        end
    end
    return breakable
end

function count_bricks(file)
    bricks = read_bricks(file)
    sort!(bricks, by=b->(b.r1[3], b.r0[3]))
    drop_bricks!(bricks)
    count(compute_breakable(bricks))
end
@show count_bricks(ARGS[1])

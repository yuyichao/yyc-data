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

mutable struct BreakSupport
    nsupport::Int
    ndropped::Int
    supports::Vector{BreakSupport}
    BreakSupport() = new(0, 0, BreakSupport[])
end

function compute_supports(bricks)
    nbricks = length(bricks)
    supports = [BreakSupport() for i in 1:nbricks]
    for i in 1:nbricks
        brick = bricks[i]
        support = supports[i]
        for j in i - 1:-1:1
            brick2 = bricks[j]
            if brick2.r1[3] >= brick.r0[3]
                continue
            elseif brick2.r1[3] < brick.r0[3] - 1
                break
            end
            @assert brick2.r1[3] == brick.r0[3] - 1
            if xy_overlap(brick, brick2)
                support2 = supports[j]
                support.nsupport += 1
                push!(support2.supports, support)
            end
        end
    end
    return supports
end

function count_fail(supports, i)
    for support in supports
        support.ndropped = 0
    end
    return count_fail(supports[i]) - 1
end

function count_fail(support::BreakSupport)
    c = 1
    for support2 in support.supports
        support2.ndropped += 1
        @assert support2.ndropped <= support2.nsupport
        if support2.ndropped == support2.nsupport
            c += count_fail(support2)
        end
    end
    return c
end

function count_fail(file)
    bricks = read_bricks(file)
    sort!(bricks, by=b->(b.r1[3], b.r0[3]))
    drop_bricks!(bricks)
    supports = compute_supports(bricks)
    return sum(count_fail(supports, i) for i in 1:length(bricks))
end
@show count_fail(ARGS[1])

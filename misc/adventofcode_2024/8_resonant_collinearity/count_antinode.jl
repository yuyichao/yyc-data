#!/usr/bin/julia

function count_antinode(file)
    lines = readlines(file)
    nrow = length(lines)
    ncol = length(lines[1])

    antennas = Dict{Char,Vector{NTuple{2,Int}}}()

    for row in 1:nrow
        line = lines[row]
        for col in 1:ncol
            c = line[col]
            if c == '.'
                continue
            end
            push!(get!(()->NTuple{2,Int}[], antennas, c), (row, col))
        end
    end

    antinodes = Set{NTuple{2,Int}}()
    for apos in values(antennas)
        n = length(apos)
        for j in 2:n
            p1 = apos[j]
            for i in 1:j - 1
                p2 = apos[i]
                ap1 = p2 .* 2 .- p1
                ap2 = p1 .* 2 .- p2
                if 1 <= ap1[1] <= nrow && 1 <= ap1[2] <= ncol
                    push!(antinodes, ap1)
                end
                if 1 <= ap2[1] <= nrow && 1 <= ap2[2] <= ncol
                    push!(antinodes, ap2)
                end
            end
        end
    end
    return length(antinodes)
end

@show count_antinode(ARGS[1])

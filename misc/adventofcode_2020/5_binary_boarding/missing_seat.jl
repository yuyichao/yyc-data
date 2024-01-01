#!/usr/bin/julia

function empty_seat(file)
    seats = Int[]
    for line in eachline(file)
        row = parse(Int, replace(line[1:7], "F"=>"0", "B"=>"1"), base=2)
        col = parse(Int, replace(line[8:10], "L"=>"0", "R"=>"1"), base=2)
        cid = row * 8 + col
        push!(seats, cid)
    end
    sort!(seats)
    nc = length(seats)
    for i in 2:nc
        if seats[i] - seats[i - 1] != 1
            @assert seats[i] - seats[i - 1] == 2
            return seats[i] - 1
        end
    end
end

@show empty_seat(ARGS[1])

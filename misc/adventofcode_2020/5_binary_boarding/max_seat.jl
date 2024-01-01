#!/usr/bin/julia

function max_seat(file)
    m = 0
    for line in eachline(file)
        row = parse(Int, replace(line[1:7], "F"=>"0", "B"=>"1"), base=2)
        col = parse(Int, replace(line[8:10], "L"=>"0", "R"=>"1"), base=2)
        cid = row * 8 + col
        if cid > m
            m = cid
        end
    end
    return m
end

@show max_seat(ARGS[1])

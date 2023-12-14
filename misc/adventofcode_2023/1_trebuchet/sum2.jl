#!/usr/bin/julia

function sum_file(file)
    values = Dict("0"=>0,
                  "1"=>1,
                  "2"=>2,
                  "3"=>3,
                  "4"=>4,
                  "5"=>5,
                  "6"=>6,
                  "7"=>7,
                  "8"=>8,
                  "9"=>9,
                  "zero"=>0,
                  "one"=>1,
                  "two"=>2,
                  "three"=>3,
                  "four"=>4,
                  "five"=>5,
                  "six"=>6,
                  "seven"=>7,
                  "eight"=>8,
                  "nine"=>9)
    re = r"[0-9]|one|two|three|four|five|six|seven|eight|nine"
    s = 0
    for line in eachline(file)
        matches = collect(eachmatch(re, line, overlap=true))
        d1 = values[matches[1].match]
        d2 = values[matches[end].match]
        s += d1 * 10 + d2
    end
    return s
end

@show sum_file(ARGS[1])

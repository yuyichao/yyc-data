#!/usr/bin/julia

const scores = Dict("A X"=>0 + 3,
                    "A Y"=>3 + 1,
                    "A Z"=>6 + 2,
                    "B X"=>0 + 1,
                    "B Y"=>3 + 2,
                    "B Z"=>6 + 3,
                    "C X"=>0 + 2,
                    "C Y"=>3 + 3,
                    "C Z"=>6 + 1)

function count_score(file)
    return sum(scores[line] for line in eachline(file))
end

@show count_score(ARGS[1])

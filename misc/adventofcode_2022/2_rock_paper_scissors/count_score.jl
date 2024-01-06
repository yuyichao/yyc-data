#!/usr/bin/julia

const scores = Dict("A X"=>1 + 3,
                    "A Y"=>2 + 6,
                    "A Z"=>3,
                    "B X"=>1,
                    "B Y"=>2 + 3,
                    "B Z"=>3 + 6,
                    "C X"=>1 + 6,
                    "C Y"=>2,
                    "C Z"=>3 + 3)

function count_score(file)
    return sum(scores[line] for line in eachline(file))
end

@show count_score(ARGS[1])

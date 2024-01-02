#!/usr/bin/julia

const M = [0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0]
M[7, 1] = 1
M[9, 1] = 1
M[1, 2] = 1
M[2, 3] = 1
M[3, 4] = 1
M[4, 5] = 1
M[5, 6] = 1
M[6, 7] = 1
M[7, 8] = 1
M[8, 9] = 1

function count_fish(file)
    states = zeros(Int, 9)
    for d in parse.(Int, split(read(file, String), ','))
        states[d + 1] += 1
    end
    return sum(M^80 * states)
end

@show count_fish(ARGS[1])

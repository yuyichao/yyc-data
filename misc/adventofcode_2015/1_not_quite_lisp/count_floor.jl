#!/usr/bin/julia

using DataStructures

function count_floor(file)
    a = counter(strip(read(file, String)))
    return a['('] - a[')']
end

@show count_floor(ARGS[1])

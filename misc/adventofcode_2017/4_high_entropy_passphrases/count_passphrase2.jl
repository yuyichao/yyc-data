#!/usr/bin/julia

using DataStructures

function check_passphrase(line)
    words = split(line)
    return length(words) == length(Set(counter.(words)))
end

@show count(check_passphrase, eachline(ARGS[1]))

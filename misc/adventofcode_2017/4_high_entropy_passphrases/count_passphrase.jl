#!/usr/bin/julia

function check_passphrase(line)
    words = split(line)
    return length(words) == length(Set(words))
end

@show count(check_passphrase, eachline(ARGS[1]))

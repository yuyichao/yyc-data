#!/usr/bin/julia

function find_start(file)
    s = read(file, String)
    len = length(s)
    for i in 1:len - 13
        if length(Set(s[i:i + 13])) == 14
            return i + 13
        end
    end
end

@show find_start(ARGS[1])

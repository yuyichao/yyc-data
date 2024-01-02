#!/usr/bin/julia

function find_start(file)
    s = read(file, String)
    len = length(s)
    for i in 1:len - 3
        if length(Set(s[i:i + 3])) == 4
            return i + 3
        end
    end
end

@show find_start(ARGS[1])

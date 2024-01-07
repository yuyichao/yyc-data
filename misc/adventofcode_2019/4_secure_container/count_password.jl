#!/usr/bin/julia

function check_password(num)
    s = string(num)
    prev = s[1]
    has_double = false
    for i in 2:6
        c = s[i]
        if c == prev
            has_double = true
        end
        if c < prev
            return false
        end
        prev = c
    end
    return has_double
end

@show count(check_password, 245182:790572)

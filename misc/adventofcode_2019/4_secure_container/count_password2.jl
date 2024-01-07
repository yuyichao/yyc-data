#!/usr/bin/julia

function check_password(num)
    s = string(num)
    prev = s[1]
    rep_len = 1
    has_double = false
    for i in 2:6
        c = s[i]
        if c == prev
            rep_len += 1
        elseif c < prev
            return false
        else
            if rep_len == 2
                has_double = true
            end
            rep_len = 1
        end
        prev = c
    end
    return rep_len == 2 || has_double
end

@show count(check_password, 245182:790572)

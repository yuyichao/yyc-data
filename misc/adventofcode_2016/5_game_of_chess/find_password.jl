#!/usr/bin/julia

using MD5

function find_password(id)
    suffix_num = 0
    password = UInt8[]
    while length(password) < 8
        h = md5(id * string(suffix_num))
        if h[1] == 0 && h[2] == 0 && h[3] >> 4 == 0
            push!(password, h[3] & 0xf)
        end
        suffix_num += 1
    end
    return password
end

# @show find_password("abc")
@show find_password("ffykfhsq")

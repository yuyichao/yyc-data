#!/usr/bin/julia

using MD5

function find_hash(prefix)
    suffix_num = 0
    while true
        h = md5(prefix * string(suffix_num))
        if h[1] == 0 && h[2] == 0 && h[3] >> 4 == 0
            return suffix_num
        end
        suffix_num += 1
    end
end

@show find_hash("abcdef")
@show find_hash("pqrstuv")
@show find_hash("yzbqklnj")

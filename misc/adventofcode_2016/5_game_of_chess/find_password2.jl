#!/usr/bin/julia

using MD5

function find_password(id)
    suffix_num = 0
    password = zeros(UInt8, 8)
    bit_set = falses(8)
    while !all(bit_set)
        h = md5(id * string(suffix_num))
        if h[1] == 0 && h[2] == 0 && h[3] >> 4 == 0
            bit_idx = (h[3] & 0xf) + 1
            bit = h[4] >> 4
            if bit_idx <= 8
                if !bit_set[bit_idx]
                    bit_set[bit_idx] = true
                    password[bit_idx] = bit
                    @show password, suffix_num
                end
            end
        end
        suffix_num += 1
    end
    return password
end

# @show find_password("abc")
@show find_password("ffykfhsq")

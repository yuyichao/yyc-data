#!/usr/bin/julia

const code = [2,4,
              1,1,
              7,5,
              4,0,
              0,3,
              1,6,
              5,5,
              3,0]

@inline function eval_round(va)
    vb = va & 7
    vb = vb ⊻ 1
    vc = va >> vb
    vb = vb ⊻ vc
    vb = vb ⊻ 6
    return vb & 7
end

function _find_digits(out, offset, prefix)
    if offset > length(out)
        return prefix
    end
    c = out[offset]
    for d in 0:7
        new_num = d | (prefix << 3)
        o = eval_round(new_num)
        if o == c
            r = _find_digits(out, offset + 1, new_num)
            if r >= 0
                return r
            end
        end
    end
    return -1
end

function find_digits(code)
    return _find_digits(reverse(code), 1, 0)
end
@show find_digits(code)

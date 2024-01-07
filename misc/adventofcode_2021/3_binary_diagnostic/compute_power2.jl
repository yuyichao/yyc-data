#!/usr/bin/julia

function filter_binaries(strings, bit, is_max)
    count1 = count(s->s[bit] == '1', strings)
    count0 = length(strings) - count1
    if count1 == 0 || count0 == 0
        return
    end
    use1 = (count1 >= count0) ⊻ (!is_max)
    filter!(s->s[bit] == (use1 ? '1' : '0'), strings)
    return
end

function compute_power(file)
    strings1 = readlines(file)
    strings2 = copy(strings1)

    nbits = length(first(strings1))
    for i in 1:nbits
        filter_binaries(strings1, i, false)
        filter_binaries(strings2, i, true)
    end
    @assert length(strings1) == 1
    @assert length(strings2) == 1
    return parse(Int, strings1[1], base=2) * parse(Int, strings2[1], base=2)
end

@show compute_power(ARGS[1])

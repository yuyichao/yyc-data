#!/usr/bin/julia

using Combinatorics

const original = Dict(Set("abcefg")=>'0',
                      Set("cf")=>'1',
                      Set("acdeg")=>'2',
                      Set("acdfg")=>'3',
                      Set("bcdf")=>'4',
                      Set("abdfg")=>'5',
                      Set("abdefg")=>'6',
                      Set("acf")=>'7',
                      Set("abcdefg")=>'8',
                      Set("abcdfg")=>'9')

function generate_permutations()
    all_perm = Dict{Set{Set{Char}},Dict{Set{Char},Char}}()
    for perm in permutations("abcdefg")
        cmap = Dict(zip("abcdefg", perm))
        keys = Set{Set{Char}}()
        val_map = Dict{Set{Char},Char}()
        for (orig_key, val) in original
            key = Set(cmap[c] for c in orig_key)
            push!(keys, key)
            val_map[key] = val
        end
        all_perm[keys] = val_map
    end
    return all_perm
end
const all_perm = generate_permutations()

function decode(file)
    c = 0
    for line in eachline(file)
        all_digits, info_digits = split(line, '|')
        digit_key = Set(Set(digit) for digit in split(all_digits))
        val_map = all_perm[digit_key]

        c += parse(Int, String(Char[val_map[Set(digit)]
                                    for digit in split(info_digits)]))
    end
    return c
end

@show decode(ARGS[1])

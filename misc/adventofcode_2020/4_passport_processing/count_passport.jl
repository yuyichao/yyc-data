#!/usr/bin/julia

function valid_passport(fields)
    for field in ("byr", "iyr", "eyr", "hgt", "hcl", "ecl", "pid")
        if !(field in fields)
            return false
        end
    end
    return true
end

function count_passport(file)
    c = 0
    fields = Set{String}()
    for line in eachline(file)
        if isempty(line)
            c += valid_passport(fields)
            empty!(fields)
            continue
        end
        for entry in split(line)
            push!(fields, match(r"(.+):", entry)[1])
        end
    end
    c += valid_passport(fields)
    return c
end

@show count_passport(ARGS[1])

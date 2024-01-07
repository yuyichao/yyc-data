#!/usr/bin/julia

function valid_passport(fields)
    for field in ("byr", "iyr", "eyr", "hgt", "hcl", "ecl", "pid")
        if !haskey(fields, field)
            return false
        end
    end
    byr = fields["byr"]
    if match(r"^\d\d\d\d$", byr) === nothing || !(1920 <= parse(Int, byr) <= 2002)
        return false
    end
    iyr = fields["iyr"]
    if match(r"^\d\d\d\d$", iyr) === nothing || !(2010 <= parse(Int, iyr) <= 2020)
        return false
    end
    eyr = fields["eyr"]
    if match(r"^\d\d\d\d$", eyr) === nothing || !(2020 <= parse(Int, eyr) <= 2030)
        return false
    end
    hgt = fields["hgt"]
    m = match(r"^(\d+)(in|cm)$", hgt)
    if m === nothing
        return false
    end
    hgt = parse(Int, m[1])
    if m[2] == "in" && !(59 <= hgt <= 76)
        return false
    end
    if m[2] == "cm" && !(150 <= hgt <= 193)
        return false
    end
    hcl = fields["hcl"]
    if match(r"^#[0-9a-f]{6}$", hcl) === nothing
        return false
    end
    ecl = fields["ecl"]
    if match(r"^(amb|blu|brn|gry|grn|hzl|oth)$", ecl) === nothing
        return false
    end
    pid = fields["pid"]
    if match(r"^\d{9}$", pid) === nothing
        return false
    end
    return true
end

function count_passport(file)
    c = 0
    fields = Dict{String,String}()
    for line in eachline(file)
        if isempty(line)
            c += valid_passport(fields)
            empty!(fields)
            continue
        end
        for entry in split(line)
            m = match(r"(.+):(.+)", entry)
            fields[m[1]] = m[2]
        end
    end
    c += valid_passport(fields)
    return c
end

@show count_passport(ARGS[1])

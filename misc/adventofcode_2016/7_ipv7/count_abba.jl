#!/usr/bin/julia

function contains_abba(str)
    for m in eachmatch(r"([a-z])([a-z])\2\1", str, overlap=true)
        if m[1] != m[2]
            return true
        end
    end
    return false
end

function supports_tls(str)
    for m in eachmatch(r"\[([^][]+)\]", str, overlap=true)
        if contains_abba(m[1])
            return false
        end
    end
    return contains_abba(str)
end

@show count(supports_tls, eachline(ARGS[1]))

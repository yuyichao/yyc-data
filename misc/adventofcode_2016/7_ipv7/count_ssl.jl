#!/usr/bin/julia

function find_aba(set, str, rev=false)
    for m in eachmatch(r"([a-z])([a-z])\1", str, overlap=true)
        if m[1] != m[2]
            if rev
                push!(set, (m[2][1], m[1][1]))
            else
                push!(set, (m[1][1], m[2][1]))
            end
        end
    end
end

function supports_ssl(str)
    segs = split(str, '[')
    out_ab = Set{NTuple{2,Char}}()
    in_ab = Set{NTuple{2,Char}}()

    find_aba(out_ab, popfirst!(segs))

    for seg in segs
        s1, s2 = split(seg, ']')
        find_aba(in_ab, s1, true)
        find_aba(out_ab, s2)
    end

    return !isempty(intersect(in_ab, out_ab))
end

@show count(supports_ssl, eachline(ARGS[1]))

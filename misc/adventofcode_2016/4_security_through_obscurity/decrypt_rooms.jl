#!/usr/bin/julia

using DataStructures

function room_name(line)
    m = match(r"(.*)-(\d+)\[(.....)\]", line)
    a = counter(m[1])
    reset!(a, '-')
    counts = collect(a)
    sort!(counts, by=x->(-x[2], x[1]))
    if String([counts[i][1] for i in 1:5]) != m[3]
        return nothing, 0
    end
    id = parse(Int, m[2])
    return String([begin
                       if c == '-'
                           ' '
                       else
                           'a' + (c - 'a' + id) % 26
                       end
                   end for c in m[1]]), id
end

function show_room_names(file)
    for line in eachline(file)
        name, id = room_name(line)
        if name !== nothing
            println((name, id))
        end
    end
end

@show show_room_names(ARGS[1])

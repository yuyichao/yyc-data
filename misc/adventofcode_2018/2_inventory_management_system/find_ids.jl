#!/usr/bin/julia

function check_ids(id1, id2)
    if id1 == id2
        return
    end
    diff_idx = 0
    for (i, (c1, c2)) in enumerate(zip(id1, id2))
        if c1 != c2
            if diff_idx != 0
                return
            end
            diff_idx = i
        end
    end
    return id1[1:diff_idx - 1] * id1[diff_idx + 1:end]
end

function find_ids(file)
    ids = readlines(file)
    for id1 in ids
        for id2 in ids
            id0 = check_ids(id1, id2)
            if id0 !== nothing
                return id0
            end
        end
    end
end

@show find_ids(ARGS[1])

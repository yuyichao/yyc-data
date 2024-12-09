#!/usr/bin/julia

function calc_checksum(file)
    s = read(file, String)
    f = Int[]
    is_blank = false
    id = 0
    for c in s
        l = c - '0'
        c = is_blank ? -1 : id
        if is_blank
            c = -1
            is_blank = false
        else
            c = id
            id += 1
            is_blank = true
        end
        for _ in 1:l
            push!(f, c)
        end
    end

    blk = 1
    function move(i)
        c = f[i]
        while blk < i
            if f[blk] == -1
                f[blk] = c
                f[i] = -1
                blk += 1
                return true
            end
            blk += 1
        end
        return false
    end

    for i in length(f):-1:1
        if f[i] == -1
            continue
        end
        if !move(i)
            break
        end
    end

    s = 0
    for i in 1:length(f)
        if f[i] == -1
            break
        end
        s += (i - 1) * f[i]
    end
    return s
end

@show calc_checksum(ARGS[1])

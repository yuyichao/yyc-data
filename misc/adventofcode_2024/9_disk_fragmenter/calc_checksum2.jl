#!/usr/bin/julia

function calc_checksum(file)
    s = read(file, String)
    f = Int[]
    is_blank = false
    id = 0
    file_info = NTuple{2,Int}[]
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
            push!(file_info, (length(f) + 1, l))
        end
        for _ in 1:l
            push!(f, c)
        end
    end

    function move_file(id)
        pos, len = file_info[id]
        if len == 0
            return
        end
        gap_pos = 0
        gap_len = 0
        for i in 1:pos - 1
            if f[i] != -1
                gap_len = 0
                continue
            end
            if gap_len == 0
                gap_pos = i
            end
            gap_len += 1
            if gap_len >= len
                for j in 1:len
                    f[pos + j - 1] = -1
                    f[gap_pos + j - 1] = id - 1
                end
                return
            end
        end
    end

    for id in length(file_info):-1:1
        move_file(id)
    end

    s = 0
    for i in 1:length(f)
        if f[i] == -1
            continue
        end
        s += (i - 1) * f[i]
    end
    return s
end

@show calc_checksum(ARGS[1])

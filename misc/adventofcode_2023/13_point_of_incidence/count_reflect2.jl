#!/usr/bin/julia

function check_col_reflect(M, col)
    ncol = size(M, 2)
    nrow = size(M, 1)
    diff_count = 0
    for i1 in 1:col
        i2 = 2 * col - i1 + 1
        if i2 > ncol
            continue
        end
        for j in 1:nrow
            if M[j, i1] != M[j, i2]
                diff_count += 1
                if diff_count >= 2
                    return false
                end
            end
        end
    end
    return diff_count == 1
end

function find_col_reflect(M)
    ncol = size(M, 2)
    for i in 1:ncol - 1
        if check_col_reflect(M, i)
            return i
        end
    end
    return
end

function find_reflect(M)
    col = find_col_reflect(M)
    if col === nothing
        return 100 * find_col_reflect(M')
    end
    return col
end

function process_lines(lines)
    M = Matrix{Bool}(undef, length(lines), length(lines[1]))
    for i in 1:length(lines)
        line = lines[i]
        for j in 1:length(line)
            M[i, j] = line[j] == '#'
        end
    end
    return find_reflect(M)
end

function count_reflect(file)
    s = 0
    lines_buff = String[]
    for line in eachline(file)
        if isempty(line)
            if !isempty(lines_buff)
                s += process_lines(lines_buff)
                empty!(lines_buff)
            end
        else
            push!(lines_buff, line)
        end
    end
    if !isempty(lines_buff)
        s += process_lines(lines_buff)
        empty!(lines_buff)
    end
    return s
end

@show count_reflect(ARGS[1])

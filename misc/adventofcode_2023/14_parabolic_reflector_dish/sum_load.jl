#!/usr/bin/julia

function load_lines(lines)
    M = Matrix{Int8}(undef, length(lines), length(lines[1]))
    for i in 1:length(lines)
        line = lines[i]
        for j in 1:length(line)
            c = line[j]
            M[i, j] = c == 'O' ? 1 : (c == '#' ? -1 : 0)
        end
    end
    return M
end

@inline function slide_north!(M)
    @inbounds for j in 1:size(M, 2)
        last_empty = 1
        for i in 1:size(M, 1)
            c = M[i, j]
            if c == 0
                continue
            elseif c == 1
                if i != last_empty
                    M[i, j] = 0
                    M[last_empty, j] = c
                end
                last_empty += 1
            else
                last_empty = i + 1
            end
        end
    end
end
# function slide_north!(M)
#     for j in 1:size(M, 2)
#         last_empty = 1
#         for i in 1:size(M, 1)
#             c = M[i, j]
#             if c == 0
#                 continue
#             elseif c == 1
#                 M[last_empty, j] = c
#                 last_empty += 1
#             else
#                 @assert c == -1
#                 M[last_empty:i - 1, j] .= 0
#                 last_empty = i + 1
#             end
#         end
#         M[last_empty:end, j] .= 0
#     end
# end

function compute_load(M)
    load = 0
    for j in 1:size(M, 2)
        for i in 1:size(M, 1)
            if M[i, j] == 1
                load += size(M, 1) - i + 1
            end
        end
    end
    return load
end

function sum_load(file)
    M = load_lines(readlines(file))
    slide_north!(M)
    return compute_load(M)
end

@show sum_load(ARGS[1])

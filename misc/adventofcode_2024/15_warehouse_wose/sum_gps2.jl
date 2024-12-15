#!/usr/bin/julia

function sum_gps(input, input_map)
    lines = readlines(input_map)
    nrow = length(lines)
    ncol = length(lines[1])
    M = Matrix{Int}(undef, nrow, ncol * 2)
    p = (0, 0)
    for row in 1:nrow
        line = lines[row]
        for col in 1:ncol
            c = line[col]
            if c == '#'
                M[row, col * 2 - 1] = -1
                M[row, col * 2] = -1
            elseif c == 'O'
                M[row, col * 2 - 1] = 1
                M[row, col * 2] = 0
            elseif c == '.'
                M[row, col * 2 - 1] = 0
                M[row, col * 2] = 0
            elseif c == '@'
                M[row, col * 2 - 1] = 0
                M[row, col * 2] = 0
                p = (row, col * 2 - 1)
            else
                error("Unknown character $c")
            end
        end
    end

    function check_vert(p, dir)
        if M[p[1], p[2]] == -1
            return false
        end
        if M[p[1], p[2]] == 1
            return check_vert(p .+ dir, dir) && check_vert(p .+ dir .+ (0, 1), dir)
        end
        @assert M[p[1], p[2]] == 0
        if M[p[1], p[2] - 1] == 1
            return check_vert(p .+ (0, -1) .+ dir, dir) && check_vert(p .+ dir, dir)
        end
        return true
    end

    function push_vert(p, dir)
        @assert M[p[1], p[2]] != -1
        if M[p[1], p[2]] == 1
            push_vert(p .+ dir, dir)
            push_vert(p .+ dir .+ (0, 1), dir)
            M[p...] = 0
            M[(p .+ dir)...] = 1
            return
        end
        @assert M[p[1], p[2]] == 0
        if M[p[1], p[2] - 1] == 1
            push_vert(p .+ (0, -1) .+ dir, dir)
            push_vert(p .+ dir, dir)
            M[p[1], p[2] - 1] = 0
            M[p[1] + dir[1], p[2] + dir[2] - 1] = 1
            return
        end
    end

    for line in eachline(input)
        for c in line
            if c == '^'
                dir = (-1, 0)
                newp = p .+ dir
                if !check_vert(newp, dir)
                    continue
                end
                push_vert(newp, dir)
                p = newp
            elseif c == 'v'
                dir = (1, 0)
                newp = p .+ dir
                if !check_vert(newp, dir)
                    continue
                end
                push_vert(newp, dir)
                p = newp
            elseif c == '>'
                dir = (0, 1)
                newp = p .+ dir
                n = 0
                oldline = M[p[1], :]
                @assert oldline[p[2]] == 0
                oldline[p[2]] = 2
                while M[(p .+ dir .* (n * 2 + 1))...] == 1
                    @assert M[(p .+ dir .* (n * 2 + 2))...] == 0
                    n += 1
                end
                if M[(p .+ dir .* (n * 2 + 1))...] == -1
                    # print("Old: ")
                    # for v in oldline
                    #     if v == -1
                    #         print('#')
                    #     elseif v == 0
                    #         print('.')
                    #     elseif v == 1
                    #         print('O')
                    #     elseif v == 2
                    #         print('@')
                    #     else
                    #         error("...")
                    #     end
                    # end
                    # println()
                    continue
                end
                @assert M[(p .+ dir .* (n * 2 + 1))...] == 0
                if n > 0
                    M[(p .+ dir)...] = 0
                    for i in 1:n
                        M[(p .+ dir .* (i * 2 - 1))...] = 0
                        M[(p .+ dir .* (i * 2))...] = 1
                    end
                end
                newline = M[p[1], :]
                @assert newline[newp[2]] == 0
                newline[newp[2]] = 2

                # print("Old: ")
                # for v in oldline
                #     if v == -1
                #         print('#')
                #     elseif v == 0
                #         print('.')
                #     elseif v == 1
                #         print('O')
                #     elseif v == 2
                #         print('@')
                #     else
                #         error("...")
                #     end
                # end
                # println()
                # print("New: ")
                # for v in newline
                #     if v == -1
                #         print('#')
                #     elseif v == 0
                #         print('.')
                #     elseif v == 1
                #         print('O')
                #     elseif v == 2
                #         print('@')
                #     else
                #         error("...")
                #     end
                # end
                # println()
                p = newp
            elseif c == '<'
                dir = (0, -1)
                newp = p .+ dir
                n = 0
                oldline = M[p[1], :]
                @assert oldline[p[2]] == 0
                oldline[p[2]] = 2
                while M[(p .+ dir .* (n * 2 + 2))...] == 1
                    @assert M[(p .+ dir .* (n * 2 + 1))...] == 0
                    n += 1
                end
                if M[(p .+ dir .* (n * 2 + 1))...] == -1
                    # print("Old: ")
                    # for v in oldline
                    #     if v == -1
                    #         print('#')
                    #     elseif v == 0
                    #         print('.')
                    #     elseif v == 1
                    #         print('O')
                    #     elseif v == 2
                    #         print('@')
                    #     else
                    #         error("...")
                    #     end
                    # end
                    # println()
                    continue
                end
                @assert M[(p .+ dir .* (n * 2 + 1))...] == 0
                if n > 0
                    @assert M[(p .+ dir)...] == 0
                    for i in 1:n
                        M[(p .+ dir .* (i * 2))...] = 0
                        M[(p .+ dir .* (i * 2 + 1))...] = 1
                    end
                end
                newline = M[p[1], :]
                @assert newline[newp[2]] == 0
                newline[newp[2]] = 2

                # print("Old: ")
                # for v in oldline
                #     if v == -1
                #         print('#')
                #     elseif v == 0
                #         print('.')
                #     elseif v == 1
                #         print('O')
                #     elseif v == 2
                #         print('@')
                #     else
                #         error("...")
                #     end
                # end
                # println()
                # print("New: ")
                # for v in newline
                #     if v == -1
                #         print('#')
                #     elseif v == 0
                #         print('.')
                #     elseif v == 1
                #         print('O')
                #     elseif v == 2
                #         print('@')
                #     else
                #         error("...")
                #     end
                # end
                # println()
                p = newp
            else
                error("Unknown character $c")
            end
        end
    end

    s = 0
    for row in 1:nrow
        for col in 1:ncol * 2
            if M[row, col] == 1
                s += (row - 1) * 100 + (col - 1)
            end
        end
    end
    return s
end

@show sum_gps(ARGS[1], ARGS[2])

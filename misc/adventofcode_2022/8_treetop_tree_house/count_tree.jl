#!/usr/bin/julia

function count_tree(file)
    lines = readlines(file)

    nrows = length(lines)
    ncols = length(lines[1])

    M = Matrix{Int}(undef, ncols, nrows)

    for i in 1:nrows
        line = lines[i]
        for j in 1:ncols
            M[i, j] = line[j] - '0'
        end
    end

    visible = zeros(Bool, ncols, nrows)

    for i in 1:nrows
        threshold = -1
        for j in 1:ncols
            h = M[i, j]
            if h > threshold
                visible[i, j] = true
                threshold = h
                if h == 9
                    break
                end
            end
        end

        threshold = -1
        for j in ncols:-1:1
            h = M[i, j]
            if h > threshold
                visible[i, j] = true
                threshold = h
                if h == 9
                    break
                end
            end
        end
    end

    for j in 1:ncols
        threshold = -1
        for i in 1:nrows
            h = M[i, j]
            if h > threshold
                visible[i, j] = true
                threshold = h
                if h == 9
                    break
                end
            end
        end

        threshold = -1
        for i in nrows:-1:1
            h = M[i, j]
            if h > threshold
                visible[i, j] = true
                threshold = h
                if h == 9
                    break
                end
            end
        end
    end

    return count(visible)
end

@show count_tree(ARGS[1])

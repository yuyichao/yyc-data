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

@inline function slide_south!(M)
    @inbounds for j in 1:size(M, 2)
        last_empty = size(M, 1)
        for i in size(M, 1):-1:1
            c = M[i, j]
            if c == 0
                continue
            elseif c == 1
                if i != last_empty
                    M[i, j] = 0
                    M[last_empty, j] = c
                end
                last_empty -= 1
            else
                last_empty = i - 1
            end
        end
    end
end

function slide_cycle!(M)
    slide_north!(M)
    slide_north!(M')
    slide_south!(M)
    slide_south!(M')
end

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

function record_cycles(M)
    d = Dict(M=>0)
    @time for i in 1:1000_000
        M = copy(M)
        slide_cycle!(M)
        if haskey(d, M)
            return d[M], i, M
        end
        d[M] = i
    end
    return
end

function sum_load(file)
    M = load_lines(readlines(file))
    (i1, i2, M) = record_cycles(M)
    nround = 1000_000_000
    nround = (nround - i2) % (i2 - i1)
    for i in 1:nround
        slide_cycle!(M)
    end
    return compute_load(M)
end

@show sum_load(ARGS[1])

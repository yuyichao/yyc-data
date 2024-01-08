#!/usr/bin/julia

function find_val(n)
    grid = zeros(Int, 1024, 1024)
    x0, y0 = 512, 512
    function fill_coord(x, y)
        j0, i0 = x + x0, y + y0
        s = 0
        for j in j0 - 1:j0 + 1
            for i in i0 - 1:i0 + 1
                s += grid[i, j]
            end
        end
        grid[i0, j0] = s
        if s >= n
            return s
        end
    end
    grid[y0, x0] = 1
    for ring in 1:512
        for n in 1:ring * 2
            v = fill_coord(ring, n - ring)
            if v !== nothing
                return v
            end
        end
        for n in 1:ring * 2
            v = fill_coord(ring - n, ring)
            if v !== nothing
                return v
            end
        end
        for n in 1:ring * 2
            v = fill_coord(-ring, ring - n)
            if v !== nothing
                return v
            end
        end
        for n in 1:ring * 2
            v = fill_coord(n - ring, -ring)
            if v !== nothing
                return v
            end
        end
    end
end

@show find_val(325489)

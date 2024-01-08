#!/usr/bin/julia

function move_button(x, y, inst)
    if inst == 'U' && y > 1
        y -= 1
    elseif inst == 'D' && y < 3
        y += 1
    elseif inst == 'L' && x > 1
        x -= 1
    elseif inst == 'R' && x < 3
        x += 1
    end
    return x, y
end

button_num(x, y) = LinearIndices((3, 3))[x, y]

function find_password(file)
    x, y = 2, 2
    pw = Int[]
    for line in eachline(file)
        for inst in line
            x, y = move_button(x, y, inst)
        end
        push!(pw, button_num(x, y))
    end
    return pw
end

@show find_password(ARGS[1])

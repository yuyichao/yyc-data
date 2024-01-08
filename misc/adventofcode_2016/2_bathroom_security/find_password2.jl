#!/usr/bin/julia

const pad = Char[ 0   0  '1'  0   0
                  0  '2' '3' '4'  0
                 '5' '6' '7' '8' '9'
                  0  'A' 'B' 'C'  0
                  0   0  'D'  0   0]

function move_button(x0, y0, inst)
    x, y = x0, y0
    if inst == 'U'
        y -= 1
    elseif inst == 'D'
        y += 1
    elseif inst == 'L'
        x -= 1
    elseif inst == 'R'
        x += 1
    else
        error()
    end
    if !checkbounds(Bool, pad, y, x) || pad[y, x] == '\0'
        return (x0, y0)
    end
    return x, y
end

button_num(x, y) = pad[y, x]

function find_password(file)
    x, y = 1, 3
    pw = Char[]
    for line in eachline(file)
        for inst in line
            x, y = move_button(x, y, inst)
        end
        push!(pw, button_num(x, y))
    end
    return String(pw)
end

@show find_password(ARGS[1])

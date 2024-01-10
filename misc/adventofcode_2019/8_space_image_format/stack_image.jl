#!/usr/bin/julia

function stack_image(file)
    str = strip(read(file, String))
    len = length(str)
    image = fill('2', 25, 6)
    npixels = 25 * 6
    for i in 1:(25 * 6):len
        for j in 1:npixels
            c0 = image[j]
            if c0 != '2'
                continue
            end
            c = str[i + j - 1]
            image[j] = c == '0' ? ' ' : c
        end
    end
    for i in 1:6
        println(String(image[:, i]))
    end
end

stack_image(ARGS[1])

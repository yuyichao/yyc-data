#!/usr/bin/julia

mutable struct Cubes
    blue::Int
    red::Int
    green::Int
end

function check_item(item, cube::Cubes)
    num, name = split(strip(item), " ")
    num = parse(Int, num)
    if name == "blue"
        cube.blue = max(cube.blue, num)
    elseif name == "red"
        cube.red = max(cube.red, num)
    elseif name == "green"
        cube.green = max(cube.green, num)
    else
        error("unknown color $name")
    end
end

function check_result(result, cube)
    for s in split(result, ",")
        check_item(s, cube)
    end
end

function count_power_file(file)
    s = 0
    for line in eachline(file)
        game_id, results = split(line, ":")
        cube = Cubes(0, 0, 0)
        for s in split(results, ";")
            check_result(s, cube)
        end
        s += cube.blue * cube.red * cube.green
    end
    return s
end

@show count_power_file(ARGS[1])

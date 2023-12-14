#!/usr/bin/julia

function check_item(item, max_blue, max_red, max_green)
    num, name = split(strip(item), " ")
    num = parse(Int, num)
    if name == "blue"
        return num <= max_blue
    elseif name == "red"
        return num <= max_red
    elseif name == "green"
        return num <= max_green
    else
        error("unknown color $name")
    end
end

function check_result(result, max_blue, max_red, max_green)
    return all(check_item(s, max_blue, max_red, max_green)
               for s in split(result, ","))
end

function count_game_file(file, max_blue, max_red, max_green)
    s = 0
    for line in eachline(file)
        game_id, results = split(line, ":")
        @assert game_id[1:4] == "Game"
        game_id = parse(Int, game_id[5:end])
        if all(check_result(s, max_blue, max_red, max_green)
               for s in split(results, ";"))
            s += game_id
        end
    end
    return s
end

@show count_game_file(ARGS[1], 14, 12, 13)

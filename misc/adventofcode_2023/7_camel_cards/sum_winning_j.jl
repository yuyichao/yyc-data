#!/usr/bin/julia

function hand_kind(hand)
    count = Dict{Int,Int}()
    jcount = 0
    for c in hand
        if c == 12
            jcount += 1
        else
            count[c] = get(count, c, 0) + 1
        end
    end
    if isempty(count)
        @assert jcount == 5
        count = [5]
    else
        count = sort!(collect(values(count)))
        count[end] += jcount
    end
    if count == [5]
        return 1
    elseif count == [1, 4]
        return 2
    elseif count == [2, 3]
        return 3
    elseif count == [1, 1, 3]
        return 4
    elseif count == [1, 2, 2]
        return 5
    elseif count == [1, 1, 1, 2]
        return 6
    elseif count == [1, 1, 1, 1, 1]
        return 7
    else
        error("...")
    end
end

function isless_card(hand1, hand2)
    k1 = hand_kind(hand1)
    k2 = hand_kind(hand2)
    if k1 != k2
        return k1 > k2
    end
    return hand1 > hand2
end

const card_order = Dict('A'=>0, 'K'=>1, 'Q'=>2, 'T'=>3, '9'=>4, '8'=>5, '7'=>6,
                        '6'=>7, '5'=>8, '4'=>9, '3'=>10, '2'=>11, 'J'=>12)

function sum_winning(file)
    hands = Tuple{Vector{Int},Int}[]
    for line in eachline(file)
        hand, bid = split(line)
        hand = [card_order[c] for c in hand]
        bid = parse(Int, bid)
        push!(hands, (hand, bid))
    end
    sort!(hands, lt=(x, y)->isless_card(x[1], y[1]))
    return sum(i * x[2] for (i, x) in enumerate(hands))
end

@show sum_winning(ARGS[1])

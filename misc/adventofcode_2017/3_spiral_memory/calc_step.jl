#!/usr/bin/julia

function calc_step(n)
    ring = ceil(Int, sqrt(n))
    if ring % 2 == 0
        ring += 1
    end
    edge_step = n - (ring - 2)^2

    edge_len = ring - 1

    edge_pos = (edge_step - 1) % edge_len + 1

    edge_mid_pos = edge_len ÷ 2

    return edge_mid_pos + abs(edge_pos - edge_mid_pos)
end

@show calc_step(325489)

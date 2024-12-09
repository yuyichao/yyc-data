#!/usr/bin/julia

function eval_logic(v0, v1, v2, v3, v4, v5, v6, v7, v8, v9)
    m1 = v0 ⊻ ~v2
    m2 = ~v0 ⊻ v3
    m3 = ~v8 ⊻ v1
    m4 = v5 ⊻ ~v3
    m5 = v2 ⊻ v4
    m6 = ~v5 ⊻ ~v1
    m7 = v6 ⊻ ~v4
    m8 = ~v6 ⊻ ~v7
    m9 = v8 ⊻ ~v9
    m10 = v7 ⊻ v9
    return ((m3 & m7 & m10),
            (m2 & m5 & m8),
            ~(m1 & m4 & m10),
            (m6 & m5 & m9),
            ~(m4 & m3 & m9),
            (m1 & m6 & m8),
            (m2 & m9 & m4),
            ~(m8 & m9 & m10))
end

function eval_logic(v)
    eval_logic((v >> 9) & 1 != 0,
               (v >> 8) & 1 != 0,
               (v >> 7) & 1 != 0,
               (v >> 6) & 1 != 0,
               (v >> 5) & 1 != 0,
               (v >> 4) & 1 != 0,
               (v >> 3) & 1 != 0,
               (v >> 2) & 1 != 0,
               (v >> 1) & 1 != 0,
               (v >> 0) & 1 != 0)
end

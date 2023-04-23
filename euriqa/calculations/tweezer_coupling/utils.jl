#!/usr/bin/julia

function add_J₊!(M::AbstractMatrix{T}, scale=1) where T
    (s1, s2) = size(M)
    @assert s1 == s2
    J = (s1 - 1) / T(2)
    for i in 1:(s1 - 1)
        m = i - J - 1
        M[i + 1, i] += sqrt((J - m) * (J + m + 1)) * scale
    end
    return M
end

function add_J₋!(M::AbstractMatrix{T}, scale=1) where T
    (s1, s2) = size(M)
    @assert s1 == s2
    J = (s1 - 1) / T(2)
    for i in 2:s1
        m = i - J - 1
        M[i - 1, i] += sqrt((J + m) * (J - m + 1)) * scale
    end
    return M
end

function add_Jz!(M::AbstractMatrix{T}, scale=1) where T
    (s1, s2) = size(M)
    @assert s1 == s2
    J = (s1 - 1) / T(2)
    for i in 1:s1
        m = i - J - 1
        M[i, i] += m * T(scale)
    end
    return M
end

# J+ = Jx + i Jy
# J- = Jx - i Jy

# (J+ + J-) / 2 = Jx
function add_Jx!(M::AbstractMatrix{T}, scale=1) where T
    add_J₊!(M, scale / T(2))
    add_J₋!(M, scale / T(2))
end

# (J+ - J-) / 2i = Jy
function add_Jy!(M::AbstractMatrix{T}, scale=1) where T
    add_J₊!(M, scale / T(2im))
    add_J₋!(M, -scale / T(2im))
end

#!/usr/bin/julia

using LinearAlgebra

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

function add_vector_shift!(M::AbstractMatrix{T}, α, u) where T
    (s1, s2) = size(M)
    @assert s1 == s2
    J = (s1 - 1) / T(2)
    V = imag(u × conj(u)) ./ 2 .* α ./ J
    add_Jx!(M, V[1])
    add_Jy!(M, V[2])
    add_Jz!(M, V[3])
    return M
end

function add_tensor_shift!(M::AbstractMatrix{T}, α, u) where T
    (s1, s2) = size(M)
    @assert s1 == s2
    J = (s1 - 1) / T(2)

    M1 = similar(M)
    M2 = similar(M)
    M1 .= zero(T)
    M2 .= zero(T)

    add_Jx!(M1, u[1])
    add_Jy!(M1, u[2])
    add_Jz!(M1, u[3])

    add_Jx!(M2, conj(u[1]))
    add_Jy!(M2, conj(u[2]))
    add_Jz!(M2, conj(u[3]))

    M′ = 3 * (M1 * M2 + M2 * M1) - 2 * J * (J + 1) * I

    M .+= M′ .* α ./ (2 * J * (2 * J - 1))

    return M
end

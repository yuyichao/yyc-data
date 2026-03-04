#!/usr/bin/julia

using AMO
using AMO: Atomic

using LinearAlgebra

const μB = 1.39962449171
const μN = 0.7622593285e-3

const g_I = -1.093 * μN / (9/2) / μB
const g_3P1 = AMO.g_sum(1, 1, AMO.g_l, 1, AMO.g_s)

ground_matrix(B) = Atomic.hyperfine_matrix(I=9/2, J=0, Bm=B, g_I=2π * g_I * μB, g_J=0,
                                           Ahf=0, Bhf=0)
excited_matrix(B) = Atomic.hyperfine_matrix(I=9/2, J=1, Bm=B, g_I=2π * g_I * μB,
                                            g_J=2π * g_3P1 * μB,
                                            Ahf=2π * -260.083, Bhf=2π * -35.355)
couple_matrix(Ωs) = Atomic.dipole_couple_matrix(0, 1, Ωs; S=9/2)

struct Kernel
    base::Matrix{ComplexF64}
    amp::Matrix{ComplexF64}
    det::Matrix{ComplexF64}

    buff::Matrix{ComplexF64}
    buff2::Matrix{ComplexF64}

    grad_amp::Matrix{ComplexF64}
    grad_det::Matrix{ComplexF64}
end

function Kernel(pol, B)
    G = ground_matrix(B)
    E = excited_matrix(B)
    C = couple_matrix(pol)

    sz_g = size(G, 1)
    sz_e = size(E, 1)
    sz = sz_g + sz_e
    @assert size(G) == (sz_g, sz_g)
    @assert size(E) == (sz_e, sz_e)
    @assert size(C) == (sz_g, sz_e)

    base = zeros(ComplexF64, sz, sz)
    amp = zeros(ComplexF64, sz, sz)
    det = zeros(ComplexF64, sz, sz)

    base[1:sz_g, 1:sz_g] .= G
    base[sz_g + 1:sz_g + sz_e, sz_g + 1:sz_g + sz_e] .= E

    amp[1:sz_g, sz_g + 1:sz_g + sz_e] .= C
    amp[sz_g + 1:sz_g + sz_e, 1:sz_g] .= C'

    for i in sz_g + 1:sz_g + sz_e
        det[i, i] = 1
    end
    return Kernel(base, amp, det, zeros(ComplexF64, sz, sz),
                  zeros(ComplexF64, sz * 2, sz * 2))
end

function get_evolution(kern::Kernel, amp, det)
    sz = size(kern.base, 1)
    buff = kern.buff
    buff .= kern.base .+ amp .* kern.amp .+ det .* kern.det

    buff2 = kern.buff2

    buff2[1:sz, 1:sz] .= buff
    buff2[sz + 1:sz * 2, sz + 1:sz * 2] .= buff
    buff2[1:sz, sz + 1:sz * 2] .= kern.amp
    buff2 = LinearAlgebra.exp!(buff2)
    kern.grad_amp .= @view(buff2[sz + 1:sz * 2, sz + 1:sz * 2])

    buff2[1:sz, 1:sz] .= buff
    buff2[sz + 1:sz * 2, sz + 1:sz * 2] .= buff
    buff2[1:sz, sz + 1:sz * 2] .= kern.det
    buff2 = LinearAlgebra.exp!(buff2)
    kern.grad_det .= @view(buff2[sz + 1:sz * 2, sz + 1:sz * 2])

    buff .= @view(buff2[1:sz, 1:sz])

    return buff, kern.grad_amp, kern.grad_det
end

#!/usr/bin/julia

using WignerSymbols
using SparseArrays

function gen_dipole(F1, F2, Ωs, δ, g1=0, g2=0)
    n1 = Int(F1 * 2 + 1)
    n2 = Int(F2 * 2 + 1)
    @assert abs(F1 - F2) <= 1

    idx1(mF1) = Int(mF1 + F1 + 1)
    idx2(mF2) = Int(mF2 + F2 + 1) + n1

    Is = Int[]
    Js = Int[]
    Vs = ComplexF64[]

    if g1 != 0
        for mF1 in -F1:F1
            if mF1 == 0
                continue
            end
            i1 = idx1(mF1)
            push!(Is, i1)
            push!(Js, i1)
            push!(Vs, g1 * mF1)
        end
    end
    if g2 != 0 || δ != 0
        for mF2 in -F2:F2
            v2 = g2 * mF2 - δ
            if v2 == 0
                continue
            end
            i2 = idx2(mF2)
            push!(Is, i2)
            push!(Js, i2)
            push!(Vs, v2)
        end
    end

    function ele(mF1, q)
        return Ωs[q + 2] * clebschgordan(Float64, F1, mF1, 1, q, F2)
    end

    function fill_q(q)
        for mF1 in max(-F1, -F2 - q):min(F1, F2 - q)
            mF2 = mF1 + q
            v = ele(mF1, q)
            if v == 0
                continue
            end

            i1 = idx1(mF1)
            i2 = idx2(mF2)
            push!(Is, i1)
            push!(Js, i2)
            push!(Vs, v)

            push!(Is, i2)
            push!(Js, i1)
            push!(Vs, conj(v))
        end
    end
    for q in -1:1
        fill_q(q)
    end
    return sparse(Is, Js, Vs, n1 + n2, n1 + n2)
end

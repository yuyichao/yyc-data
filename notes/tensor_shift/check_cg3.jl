#!/usr/bin/julia

using WignerSymbols
using RationalRoots
using HalfIntegers

function f0_kernel(j0, m0, j1, m1, j2, m2, j1′, m1′, j2′, m2′, k, q)
    return clebschgordan(j1, m1, j0, m0, j1′, m1′) * clebschgordan(j2, m2, j0, m0, j2′, m2′) * clebschgordan(j2, m2, k, q, j1, m1)
end

function f0(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
    return sum(f0_kernel(j0, m0, j1, m1, j2, m2, j1′, m1′, j2′, m2′, k, q)
               for m0 in -j0:j0, m1 in -j1:j1, m2 in -j2:j2; init=BigFloat(0))
end

function f1_kernel(j0, m0, j1, m1, j2, m2, j1′, m1′, j2′, m2′, k, q)
    return (-1)^(2 * j0 - j1 - 2 * j2 + k - m1′ - m2′ - m1) *
        sqrt(RationalRoot{BigInt}((2 * j1′ + 1) * (2 * j2′ + 1) * (2 * j1 + 1))) *
        wigner3j(j1, j0, j1′, m1, m0, -m1′) * wigner3j(j2, j0, j2′, m2, m0, -m2′) * wigner3j(j2, k, j1, m2, q, -m1)
end

function f1(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
    return sum(f1_kernel(j0, m0, j1, m1, j2, m2, j1′, m1′, j2′, m2′, k, q)
               for m0 in -j0:j0, m1 in -j1:j1, m2 in -j2:j2; init=0)
end

function f2_kernel(j0, m0, j1, m1, j2, m2, j1′, m1′, j2′, m2′, k, q)
    return (-1)^(2 * j0 - j2 + 2 * k - 2 * m2′ - 3 * m1′ + m0 + m1 + m2) *
        sqrt(RationalRoot{BigInt}((2 * j1′ + 1) * (2 * j2′ + 1) * (2 * j1 + 1))) *
        wigner3j(j1′, j1, j0, -m1′, m1, -m0) * wigner3j(j2, j2′, j0, -m2, m2′, m0) * wigner3j(j2, j1, k, m2, -m1, q)
end

function f2(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
    return sum(f2_kernel(j0, m0, j1, m1, j2, m2, j1′, m1′, j2′, m2′, k, q)
               for m0 in -j0:j0, m1 in -j1:j1, m2 in -j2:j2; init=0)
end

function f3(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
    return (-1)^(j0 - j1 - 2 * j2 + 2 * k - m2′ + q) *
        sqrt(RationalRoot{BigInt}((2 * j1′ + 1) * (2 * j2′ + 1) * (2 * j1 + 1))) *
        wigner3j(j1′, j2′, k, -m1′, m2′, q) * wigner6j(j1, j2, k, j2′, j1′, j0)
end

function f4(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
    return (-1)^(j0 + j1 + j2′ + k) *
        sqrt(RationalRoot{BigInt}((2 * j1 + 1) * (2 * j2′ + 1))) *
        clebschgordan(j2′, m2′, k, q, j1′, m1′) * wigner6j(j1, j2, k, j2′, j1′, j0)
end

function check(j0, j1, j2, k)
    for j1′ in abs(j1 - j0):(j1 + j0)
        for j2′ in abs(j2 - j0):(j2 + j0)
            for m1′ in -j1′:j1′
                for m2′ in -j2′:j2′
                    for q in -k:k
                        v0 = f0(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
                        # v1 = f1(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
                        # v2 = f2(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
                        # v3 = f3(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
                        v4 = f4(j0, j1, j2, j1′, m1′, j2′, m2′, k, q)
                        # if abs(v0 - v1) > 1e-10
                        #     @show j0, j1, j2, j1′, m1′, j2′, m2′, k, q
                        #     @show v0, v1, v0 - v1
                        # end
                        # if abs(v0 - v2) > 1e-10
                        #     @show j0, j1, j2, j1′, m1′, j2′, m2′, k, q
                        #     @show v0, v2, v0 - v2
                        # end
                        # if abs(v0 - v3) > 1e-10
                        #     @show j0, j1, j2, j1′, m1′, j2′, m2′, k, q
                        #     @show v0, v3, v0 - v3
                        # end
                        if abs(v0 - v4) > 1e-10
                            @show j0, j1, j2, j1′, m1′, j2′, m2′, k, q
                            @show v0, v4, v0 - v4
                        end
                    end
                end
            end
        end
    end
end

function check(jmax, kmax)
    for dk in 0:kmax * 2
        k = half(dk)
        for dj0 in 0:jmax * 2
            j0 = half(dj0)
            for dj1 in 0:jmax * 2
                j1 = half(dj1)
                for dj2 in 0:jmax * 2
                    j2 = half(dj2)
                    if abs(j1 - j2) > k || (dj1 - dj2 - dk) % 2 != 0
                        continue
                    end
                    check(j0, j1, j2, k)
                end
            end
        end
    end
end

#!/usr/bin/julia

using CGcoefficient

function f0_kernel(k, p, F, mF, F′, mF′, F′′, mF′′, q, q′)
    return ((-1)^q * CG(1, 1, k, q, -q′, p) * CG(F′′, 1, F, mF′′, -q, mF) *
        CG(F′′, 1, F′, mF′′, -q′, mF′))
end

function f0(k, p, F, mF, F′, mF′, F′′)
    return sum(f0_kernel(k, p, F, mF, F′, mF′, F′′, mF′′, q, q′)
               for mF′′ in (-F′′:F′′), q in -1:1, q′ in -1:1;
                   init=exact_sqrt(0))
end

function f1_kernel(k, p, F, mF, F′, mF′, F′′, mF′′, q, q′)
    return ((-1)^q * threeJ(1, 1, k, q, -q′, -p) * threeJ(F′′, 1, F, mF′′, -q, -mF) *
        threeJ(F′′, 1, F′, mF′′, -q′, -mF′))
end

function f1(k, p, F, mF, F′, mF′, F′′)
    return ((-1)^Int(2 * F′′ + p - mF - mF′) *
        exact_sqrt((2 * k + 1) * (2 * F + 1) * (2 * F′ + 1)) *
        sum(f1_kernel(k, p, F, mF, F′, mF′, F′′, mF′′, q, q′)
            for mF′′ in (-F′′:F′′), q in -1:1, q′ in -1:1;
                init=exact_sqrt(0)))
end

function f2_kernel(k, p, F, mF, F′, mF′, F′′, mF′′, q, q′)
    return ((-1)^(-q) * threeJ(1, k, 1, q′, -p, -q) * threeJ(1, F, F′′, q, -mF, -mF′′) *
        threeJ(F′′, F′, 1, mF′′, mF′, -q′))
end

function f2(k, p, F, mF, F′, mF′, F′′)
    return ((-1)^Int(2 * F′′ + p - mF - mF′) *
        exact_sqrt((2 * k + 1) * (2 * F + 1) * (2 * F′ + 1)) *
        sum(f2_kernel(k, p, F, mF, F′, mF′, F′′, mF′′, q, q′)
            for mF′′ in (-F′′:F′′), q in -1:1, q′ in -1:1;
                init=exact_sqrt(0)))
end

function f_final(k, p, F, mF, F′, mF′, F′′)
    return ((-1)^Int(F′′ + mF - p) * exact_sqrt((2 * k + 1) * (2 * F + 1) * (2 * F′ + 1)) *
        threeJ(k, F, F′, p, mF, -mF′) * sixJ(k, F, F′, F′′, 1, 1))
end

function check_upto(Fmax, F′max)
    for k in 0:2
        for dF in 0:2 * Fmax
            F = dF//2
            for dF′ in (dF % 2):2:2 * F′max
                F′ = dF′//2
                if abs(F - F′) > k
                    continue
                end
                for F′′ in max(F, F′) - 1:min(F, F′) + 1
                    if abs(F - F′′) > 1 || abs(F′ - F′′) > 1
                        continue
                    end
                    for p in -k:k
                        for mF in -F:F
                            for mF′ in -F′:F′
                                if mF′ != p + mF
                                    continue
                                end
                                v1 = f0(k, p, F, mF, F′, mF′, F′′)
                                v2 = f_final(k, p, F, mF, F′, mF′, F′′)
                                if v1 - v2 != 0
                                    @assert v1 + v2 == 0
                                    @show k, p, F, mF, F′, mF′, F′′
                                    @show v1, v2, v1 - v2
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

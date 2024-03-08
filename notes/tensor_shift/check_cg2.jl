#!/usr/bin/julia

using CGcoefficient

function f0_kernel(j, p, m, m′, q, q′, m′′)
    return CG(1, 1, 2, q, q′, p) * CG(j, 1, j, m′′, q, m) * CG(j, 1, j, m′, q′, m′′)
end

function f0(j, p, m, m′)
    return sum(f0_kernel(j, p, m, m′, q, q′, m′′)
               for m′′ in -j:j, q in -1:1, q′ in -1:1;
                   init=exact_sqrt(0))
end


function f_final(j, p, m, m′)
    return (exact_sqrt(5) * (2j + 1) * (-1)^Int(j - m) *
        sixJ(2, j, j, j, 1, 1) * threeJ(j, 2, j, m′, p, -m))
end

function f_final_cg(j, p, m, m′)
    return (exact_sqrt(5 * (2j + 1)) * (-1)^Int(2j) *
        sixJ(2, j, j, j, 1, 1) * CG(j, 2, j, m′, p, m))
end

function check_upto(jmax)
    for dj in 0:2 * jmax
        j = dj//2
        for p in -2:2
            for m in -j:j
                for m′ in -j:j
                    v1 = f0(j, p, m, m′)
                    v2 = f_final(j, p, m, m′)
                    v3 = f_final_cg(j, p, m, m′)
                    if v1 - v2 != 0
                        @show j, p, m, m′
                        @show v1, v2, v1 - v2
                    end
                    if v1 - v3 != 0
                        @show j, p, m, m′
                        @show v1, v3, v1 - v3
                    end
                end
            end
        end
    end
end

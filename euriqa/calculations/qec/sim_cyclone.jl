#!/usr/bin/julia

include("sim_cyclone_lib.jl")

using NaCsCalc.Format: Unc

function parse_lists(indexes0, fidelities)
    pd2s = Float64[]
    prev_fid = 1.0
    for i in 1:length(indexes0)
        index = indexes0[i] + 1
        fid = fidelities[index]
        p = prev_fid - fid
        prev_fid = fid
        @assert p >= 0
        push!(pd2s, p)
    end
    return pd2s
end

function sim_513(pds, pd2s, n, final_round)
    stabs_x = [[true, false, false, true, true],
               [false, true, false, false, true],
               [false, false, true, true, false],
               [false, false, false, false, false]]
    stabs_z = [[true, false, true, false, true],
               [true, false, false, true, true],
               [true, false, false, true, true],
               [false, true, true, true, true]]

    logics_x = [[true, false, true, true, false],
                [false, false, false, false, false]]
    logics_z = [[false, false, false, false, false],
                [true, false, true, true, false]]

    stab_orders = [[1, 4, 5, 3],
                   [2, 5, 1, 4],
                   [3, 4, 1, 5],
                   [2, 3, 4, 5]]
    # stab_orders = nothing

    T = UInt128
    init = UniformInit{T}(0.0)
    rngs = RNGs{T}(zeros(4), pds, pd2s)
    circ = RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
                                    rngs, init; stab_orders=stab_orders)
    err_stat = run(circ, n)
    return Unc(err_stat.err / err_stat.tot,
               sqrt(err_stat.err) / err_stat.tot)
end

const base_fidelities = [0.9918, 0.9676012390492177, 0.9676012390492177,
                         0.9676012390492177, 0.9592798683933944, 0.9502626376304966,
                         0.9413301688367699, 0.9150334040229514, 0.9150334040229514,
                         0.9150334040229514, 0.9058464686465609, 0.8788353885709433,
                         0.8788353885709433, 0.8788353885709433, 0.8690557083669258,
                         0.840364806754743, 0.840364806754743, 0.840364806754743,
                         0.830400097022168, 0.8012093254309735, 0.8012093254309735,
                         0.8012093254309735, 0.7919745226493101, 0.7817787060026847,
                         0.7727678620950013, 0.7628193104598192, 0.7628193104598192,
                         0.7628193104598192, 0.7628193104598192, 0.7534167386135466,
                         0.744130063620057, 0.7169512220427975, 0.7169512220427975,
                         0.7169512220427975, 0.7084188830546323, 0.6978640537710332,
                         0.6671337758271189, 0.6671337758271189, 0.6671337758271189,
                         0.6574999029813094]
const base_indexes = sort!([0, 5, 4, 6, 10, 14, 18, 23, 22, 25, 24, 29, 30, 35, 34, 39])

const cyclone_fidelities = [0.9918, 0.98366724, 0.9756011686320001,
                            0.9756011686320001, 0.9756011686320001, 0.9756011686320001,
                            0.9756011686320001, 0.9756011686320001, 0.9756011686320001,
                            0.9756011686320001, 0.9756011686320001, 0.9756011686320001,
                            0.9756011686320001, 0.9756011686320001, 0.9756011686320001,
                            0.966820758114312, 0.9581193712912832, 0.9494962969496616,
                            0.9409508302771147, 0.9157731229582347, 0.8912691139086345,
                            0.8674207764925971, 0.8442105664262386, 0.8442105664262386,
                            0.8442105664262386, 0.8442105664262386, 0.8442105664262386,
                            0.8442105664262386, 0.8442105664262386, 0.8442105664262386,
                            0.8442105664262386, 0.8351493730132636, 0.8264081429090582,
                            0.8176482165942222, 0.8089811454983235, 0.7835272729265135,
                            0.7588742838270901, 0.7346998705895965, 0.7118709137819932,
                            0.7118709137819932, 0.7118709137819932, 0.7118709137819932,
                            0.7118709137819932, 0.7118709137819932, 0.7118709137819932,
                            0.7118709137819932, 0.7118709137819932, 0.7028380628537817,
                            0.6945601923357261, 0.6861328620020526, 0.6777620410856275,
                            0.6691619938531855]
const cyclone_indexes = sort!([0, 1, 2, 15, 16, 17, 18, 31,
                               32, 33, 34, 47, 48, 49, 50, 51])

function sim_cyclone_513()
    pds = zeros(8)
    pd2s_base = parse_lists(base_indexes, base_fidelities)
    pd2s_cyclone = parse_lists(cyclone_indexes, cyclone_fidelities)
    @show pd2s_base
    @show pd2s_cyclone
    @show pd2s_base .- pd2s_cyclone
    le_base = sim_513(pds, pd2s_base, 1000_000_000, false)
    # @show sim_513(pds, pd2s_base, 1000_000_000, true)
    le_cyclone = sim_513(pds, pd2s_cyclone, 1000_000_000, false)
    # @show sim_513(pds, pd2s_cyclone, 1000_000_000, true)
    println("Logical error:")
    println("     Base: $(le_base * 100) %")
    println("  Cyclone: $(le_cyclone * 100) %")
end

sim_cyclone_513()

#!/usr/bin/julia

include("sim_cyclone_lib.jl")

function parse_lists(indexes0, fidelitis, counts)
    pds = Float64[]
    pd2s = Float64[]

    prev_fid = 1.0
    cur_type = 1
    pos_offset = 0
    for c in counts
        for i in 1:c
            index = indexes0[i + pos_offset] + 1
            fid = fidelities[index]
            p = prev_fid - fid
            @assert p >= 0
            if cur_type == 1
                push!(pds, p)
            else
                push!(pd2s, p)
            end
            prev_fid = fid
        end
        cur_type = 3 - cur_type
        pos_offset += c
    end
    return pds, pd2s
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

    T = UInt128
    init = UniformInit{T}(0.0)
    rngs = RNGs{T}(zeros(4), pds, pd2s)
    circ = RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
                                    rngs, init, stab_orders)
    err_stat = run(circ, n)
    return (err_stat.err, sqrt(err_stat.err)) / err_stat.tot
end

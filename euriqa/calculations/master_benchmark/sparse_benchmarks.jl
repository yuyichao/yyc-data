#!/usr/bin/julia

using BenchmarkTools

const outfile = ARGS[1]

include("utils.jl")

include("qt.jl")
include("qo.jl")
include("simple.jl")
include("opt.jl")

const P1 = gen_dipole(10, 11, (2.0, 1.2im, -2.1), 0.5, 0.1, 0.05), 11, 33
const P2 = gen_dipole(10, 9, (3.2im, 2.0, -2.1), -0.5, 0.2, 0.1), 11, 12
const P3 = gen_dipole(2, 1, (3.2im, -2.1, 0.0), 0.0, 0.1, 0.1), 3, 7
const P4 = gen_dipole(9, 8, (0.0, 4.0, 0.0), 0.0), 10, 28

function test_sparse(M)
    println(M)
    println("P1")
    p = M.sparse_problem(P1...)
    t1 = @btimed($(M.psolve)($(p), $(range(0, 10, 2)))).time
    println("P2")
    p = M.sparse_problem(P2...)
    t2 = @btimed($(M.psolve)($(p), $(range(0, 10, 2)))).time
    println("P3")
    p = M.sparse_problem(P3...)
    t3 = @btimed($(M.psolve)($(p), $(range(0, 10, 2)))).time
    println("P4")
    p = M.sparse_problem(P4...)
    t4 = @btimed($(M.psolve)($(p), $(range(0, 10, 2)))).time
    return t1, t2, t3, t4
end

t_qt = test_sparse(QTTest)
t_qo = test_sparse(QOTest)
t_simple = test_sparse(SimpleTest)
t_opt = test_sparse(OPTTest)

function write_times(io, name, ts)
    println(io, name, ",", join(ts, ','))
end

open(outfile, "w") do io
    println(io, "Name,t1,t2,t3,t4")
    write_times(io, "QT", t_qt)
    write_times(io, "QO", t_qo)
    write_times(io, "Simple", t_simple)
    write_times(io, "OPT", t_opt)
end

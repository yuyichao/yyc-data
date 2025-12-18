#!/usr/bin/julia

using BenchmarkTools

const outfile = ARGS[1]

include("qt.jl")
include("qo.jl")
include("simple.jl")
include("opt.jl")

function test_motion(M)
    println(M)
    println("n=10")
    p = M.motion_problem(2π * 0.13, 2π * 0.032, 0.0, 0.0, 10)
    t1_10 = @btimed($(M.psolve)($(p), $(range(0, 4 * π / (2π * 0.13), 2)))).time
    p = M.motion2_problem(2π * 0.13, 2π * 0.032, 0.0, 0.0, 10)
    t2_10 = @btimed($(M.psolve)($(p), $(range(0, 4 * π / (2π * 0.13), 2)))).time
    println("n=50")
    p = M.motion_problem(2π * 0.13, 2π * 0.032, 0.0, 0.0, 50)
    t1_50 = @btimed($(M.psolve)($(p), $(range(0, 4 * π / (2π * 0.13), 2)))).time
    p = M.motion2_problem(2π * 0.13, 2π * 0.032, 0.0, 0.0, 50)
    t2_50 = @btimed($(M.psolve)($(p), $(range(0, 4 * π / (2π * 0.13), 2)))).time
    println("n=250")
    p = M.motion_problem(2π * 0.13, 2π * 0.032, 0.0, 0.0, 250)
    t1_250 = @btimed($(M.psolve)($(p), $(range(0, 4 * π / (2π * 0.13), 2)))).time
    p = M.motion2_problem(2π * 0.13, 2π * 0.032, 0.0, 0.0, 250)
    t2_250 = @btimed($(M.psolve)($(p), $(range(0, 4 * π / (2π * 0.13), 2)))).time
    return t1_10, t1_50, t1_250, t2_10, t2_50, t2_250
end

t_qt = test_motion(QTTest)
t_qo = test_motion(QOTest)
t_simple = test_motion(SimpleTest)
t_opt = test_motion(OPTTest)

function write_times(io, name, ts)
    println(io, name, ",", join(ts, ','))
end

open(outfile, "w") do io
    println(io, "Name,t1_10,t1_50,t1_250,t2_10,t2_50,t2_250")
    write_times(io, "QT", t_qt)
    write_times(io, "QO", t_qo)
    write_times(io, "Simple", t_simple)
    write_times(io, "OPT", t_opt)
end

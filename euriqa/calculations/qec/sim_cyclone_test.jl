#!/usr/bin/julia

include("sim_cyclone_lib.jl")

# Steane
# const stabs_x = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false]]
# const stabs_z = [[false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]

# const stabs_x = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]
# const stabs_z = [[false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]

# const stabs_x = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false]]
# const stabs_z = [[true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true],
#                  [true, true, true, true, false, false, false],
#                  [false, true, true, false, true, true, false],
#                  [false, false, true, true, false, true, true]]

# const logics_x = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, false, false, false]]
# const logics_z = [[false, false, false, false, false, false, false],
#                   [false, false, false, false, true, true, true]]

# const logics_x = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, true, true, true]]
# const logics_z = [[false, false, false, false, false, false, false],
#                   [false, false, false, false, true, true, true]]

# const logics_x = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, false, false, false]]
# const logics_z = [[false, false, false, false, true, true, true],
#                   [false, false, false, false, true, true, true]]

# Shor
# const stabs_x = [[true, true, false, false, false, false, false, false, false],
#                  [false, true, true, false, false, false, false, false, false],
#                  [false, false, false, true, true, false, false, false, false],
#                  [false, false, false, false, true, true, false, false, false],
#                  [false, false, false, false, false, false, true, true, false],
#                  [false, false, false, false, false, false, false, true, true],
#                  [false, false, false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false, false, false]]
# const stabs_z = [[false, false, false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false, false, false],
#                  [false, false, false, false, false, false, false, false, false],
#                  [true, true, true, true, true, true, false, false, false],
#                  [false, false, false, true, true, true, true, true, true]]

# const logics_x = [[true, false, false, true, false, false, true, false, false],
#                   [false, false, false, false, false, false, false, false, false]]
# const logics_z = [[false, false, false, false, false, false, false, false, false],
#                   [true, true, true, false, false, false, false, false, false]]

# 513
const stabs_x = [[true, false, false, true, true],
                 [false, true, false, false, true],
                 [false, false, true, true, false],
                 [false, false, false, false, false]]
const stabs_z = [[true, false, true, false, true],
                 [true, false, false, true, true],
                 [true, false, false, true, true],
                 [false, true, true, true, true]]

const logics_x = [[true, false, true, true, false],
                  [false, false, false, false, false]]
const logics_z = [[false, false, false, false, false],
                  [true, false, true, true, false]]

const stab_orders = [[1, 4, 5, 3],
                     [2, 5, 1, 4],
                     [3, 4, 1, 5],
                     [2, 3, 4, 5]]

function calc_errors(ps, n)
    eps = Vector{Float64}(undef, length(ps))
    T = UInt128
    nstab = length(stabs_x)
    nq = length(stabs_x[1])
    for (i, p) in enumerate(ps)
        init = UniformInit2{T}(nq, p, 0.0)
        # init = UniformInit{T}(p)
        rngs = RNGs{T}(zeros(4), zeros(4 * 2), zeros(16))
        # rngs = RNGs{T}(ones(6) .* 0.0001, ones(6 * 2) .* 0.0001, ones(24) .* 0.0001)
        # rngs = RNGs{T}(ones(6) .* p, ones(6 * 2) .* p, ones(24) .* p)
        circ = RawStabMeasureCircuit{T}(stabs_x, stabs_z, logics_x, logics_z,
                                        rngs, init; stab_orders=stab_orders)
        err_stat = run(circ, n)
        eps[i] = err_stat.err / err_stat.tot
    end
    return eps
end

const ps = range(0.0001, 1, 200)
const ps2 = exp10.(range(-5, 0, 1001))

const prefix = joinpath(@__DIR__, "imgs/steane_circuit/")

const eps2 = (
    ps = ps2,
    n0 = @time(calc_errors(ps2, 320 * 256 * 8)),
)

function plot_lines(eps)
    plot(eps.ps, eps.n0, "C1")
    plot(eps.ps, eps.ps, "C0--")
end

figure()
plot_lines(eps2)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([0, 0.15])
ylim([0, 0.15])
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)zoomin")

figure()
plot_lines(eps2)
grid()
xlabel("Physical Error")
ylabel("Logical Error")
xlim([1e-5, 0.15])
ylim([1e-7, 0.11])
xscale("log")
yscale("log")
title("[[7, 1, 3]] Circuit")
NaCsPlot.maybe_save("$(prefix)log")

NaCsPlot.maybe_show()

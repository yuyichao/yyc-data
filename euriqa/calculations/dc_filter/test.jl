#!/usr/bin/julia

function add_parallel(vi, z)
    (vout, iout) = vi
    return vout, vout / z + iout
end
function add_series(vi, z)
    (vout, iout) = vi
    return vout + iout * z, iout
end

zr(R) = R
zc(C, ω) = 1 / (im * ω * C)
zl(L, ω) = im * ω * L

function filter_435(ω, zload)
    vi = (1, 0)
    # load
    vi = add_parallel(vi, zload)

    vi = add_parallel(vi, zc(14.7e-6, ω))
    vi = add_series(vi, zr(20))
    vi = add_parallel(vi, zc(14.7e-6, ω))
    vi = add_series(vi, zr(20))
    vi = add_parallel(vi, zc(14.7e-6, ω))
    vi = add_series(vi, zr(20))
    vi = add_parallel(vi, zc(14.7e-6, ω))

    return 1 / vi[1]
end

function filter_inv_order(ω, zload)
    vi = (1, 0)
    # load
    vi = add_parallel(vi, zload)

    # stage 3.2
    vi = add_parallel(vi, zc(200e-6, ω))
    vi = add_series(vi, zl(100e-6, ω) + zr(10))
    # stage 3.1
    vi = add_parallel(vi, zc(200e-6, ω))
    vi = add_series(vi, zl(100e-6, ω) + zr(10))

    # stage 2.2
    vi = add_parallel(vi, zc(100e-6, ω))
    vi = add_series(vi, zl(10e-6, ω) + zr(10))
    # stage 2.1
    vi = add_parallel(vi, zc(100e-6, ω))
    vi = add_series(vi, zl(10e-6, ω) + zr(10))

    # stage 1
    vi = add_parallel(vi, zc(33e-9, ω))
    vi = add_series(vi, zr(150))

    return 1 / vi[1]
end

function filter(ω, zload, stage12r=10)
    vi = (1, 0)
    # load
    vi = add_parallel(vi, zload)

    # stage 3
    vi = add_parallel(vi, zc(33e-9, ω))
    vi = add_series(vi, zr(150)) # 32kHz

    # stage 2.2
    vi = add_parallel(vi, zc(100e-6, ω))
    vi = add_series(vi, zl(10e-6, ω) + zr(stage12r)) # L:5kHz, R:160Hz
    # stage 2.1
    vi = add_parallel(vi, zc(100e-6, ω))
    vi = add_series(vi, zl(10e-6, ω) + zr(stage12r))

    # stage 1.2
    vi = add_parallel(vi, zc(200e-6, ω))
    vi = add_series(vi, zl(100e-6, ω) + zr(stage12r)) # L:1.1kHz, R:80Hz
    # stage 1.1
    vi = add_parallel(vi, zc(200e-6, ω))
    vi = add_series(vi, zl(100e-6, ω) + zr(stage12r))

    return 1 / vi[1]
end

#!/usr/bin/julia

using PyPlot

function impow(i)
    i = mod(i, 4)
    if i == 0
        return complex(1, 0)
    elseif i == 1
        return complex(0, 1)
    elseif i == 2
        return complex(-1, 0)
    else
        return complex(0, -1)
    end
end

function gen_single_drives()
    base_pol = (20e6, 5e6, 20e6)
    sideband_ratio = 1 / sqrt(4)

    drives_pump = Drives()
    push!(drives_pump.d935, SingleDrive(0, base_pol))
    push!(drives_pump.d935, SingleDrive(-2π * 3.0695e9,
                                        base_pol .* (sideband_ratio * im)))
    push!(drives_pump.d935, SingleDrive(2π * 3.0695e9,
                                        base_pol .* (sideband_ratio * -im)))

    drives_detect = Drives()
    push!(drives_detect.d935, SingleDrive(0, base_pol .* sqrt(1 + 2 * sideband_ratio^2)))

    return drives_pump, drives_detect
end

function gen_multi_drives()
    base_pol = (5e6, 5e6, 5e6)
    sideband_ratio = 1 / sqrt(2)
    zeeman_det = 2π * 5.25e6

    drives_pump = Drives()
    for j in -1:1
        push!(drives_pump.d935, SingleDrive(zeeman_det * j, base_pol .* impow(j)))
        push!(drives_pump.d935, SingleDrive(-2π * 3.0695e9 + zeeman_det * j,
                                            base_pol .* (sideband_ratio * impow(j + 1))))
        push!(drives_pump.d935, SingleDrive(2π * 3.0695e9 + zeeman_det * j,
                                            base_pol .* (sideband_ratio * impow(j - 1))))
    end

    drives_detect = Drives()
    pol = base_pol .* sqrt(1 + 2 * sideband_ratio^2)
    for j in -1:1
        push!(drives_detect.d935, SingleDrive(zeeman_det * j, pol))
    end
    return drives_pump, drives_detect
end

struct PumpAndDetectResult{T1,R1,T2,R2,T3,R3,T4,R4}
    ts_detect_d20::T1
    ρs_detect_d20::R1
    ts_detect_d1::T2
    ρs_detect_d1::R2
    ts_pump_d2::T3
    ρs_pump_d2::R3
    ts_pump_d1::T4
    ρs_pump_d1::R4
end

function PumpAndDetectResult(sys, drives_pump, drives_detect)
    ρ0_d1 = (get_ρ(sys, :D, 1, -1) + get_ρ(sys, :D, 1, 0) + get_ρ(sys, :D, 1, 1)) / 3
    ρ0_d2 = (get_ρ(sys, :D, 2, -2) + get_ρ(sys, :D, 2, -1) + get_ρ(sys, :D, 2, 0) +
        get_ρ(sys, :D, 2, 1) + get_ρ(sys, :D, 2, 2)) / 5

    ts_detect_d20, ρs_detect_d20 = evolve(drives_detect, sys, get_ρ(sys, :D, 2, 0),
                                           range(0, 100e-6, 10001), maxiters=1e8)
    ts_detect_d1, ρs_detect_d1 = evolve(drives_detect, sys, ρ0_d1,
                                         range(0, 10e-6, 10001), maxiters=1e8)
    ts_pump_d2, ρs_pump_d2 = evolve(drives_pump, sys, ρ0_d2, range(0, 10e-6, 10001),
                                     maxiters=1e8)
    ts_pump_d1, ρs_pump_d1 = evolve(drives_pump, sys, ρ0_d1, range(0, 10e-6, 10001),
                                     maxiters=1e8)

    return PumpAndDetectResult(ts_detect_d20, ρs_detect_d20, ts_detect_d1, ρs_detect_d1,
                               ts_pump_d2, ρs_pump_d2, ts_pump_d1, ρs_pump_d1)
end

function plot_result(res::PumpAndDetectResult)
    figure(figsize=[6.4 * 2, 4.8 * 2])

    subplot(2, 2, 1)
    plot(res.ts_pump_d2 .* 1e6,
         [real(ρ.data[i,i]) for ρ in res.ρs_pump_d2, i in 1:20])
    grid()
    xlim([0, 5])
    ylim([0, 0.35])
    xlabel("t (\$\\mu s\$)")
    ylabel("Population")
    title("F=2 pump repump")

    subplot(2, 2, 2)
    plot(res.ts_pump_d1 .* 1e6,
         [real(ρ.data[i,i]) for ρ in res.ρs_pump_d1, i in 1:20])
    grid()
    xlim([0, 5])
    ylim([0, 0.35])
    xlabel("t (\$\\mu s\$)")
    ylabel("Population")
    title("F=1 pump repump")

    subplot(2, 2, 3)
    plot(res.ts_detect_d20 .* 1e6,
         [real(ρ.data[i,i]) for ρ in res.ρs_detect_d20, i in 1:20])
    grid()
    xlim([0, 100])
    ylim([0.998, 1])
    xlabel("t (\$\\mu s\$)")
    ylabel("Population")
    title("\$|2,0\\rangle\$ detect depump")

    subplot(2, 2, 4)
    plot(res.ts_detect_d1 .* 1e6,
         [real(ρ.data[i,i]) for ρ in res.ρs_detect_d1, i in 1:20])
    grid()
    xlim([0, 5])
    ylim([0, 0.35])
    xlabel("t (\$\\mu s\$)")
    ylabel("Population")
    title("F=1 detect repump")

    tight_layout()
end

#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
import NaCsCalc.Format: Unc, Sci
using NaCsData.Fitting: fit_data
using Images
using Statistics
using Printf

function load_sif_prefix(prefix, idxs)
    return [Float64.(load("$(prefix)_$(idx).sif").data) for idx in idxs]
end

const currents = [3, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0]
const times = [15, 8, 6, 6, 4, 4, 4, 4, 4, 4]
const freq_res = 725258845
const freq_off = 725258605

function load_images(dir, current, time)
    if isinteger(current)
        cur_str = string(round(Int, current))
    else
        cur_str = @sprintf("%.1f", current)
        cur_str = replace(cur_str, "."=>"p")
    end
    prefix = joinpath(dir, "$(time)s_$(cur_str)A_754uW")
    return (res=load_sif_prefix("$(prefix)_$(freq_res)", 0:9),
            off=load_sif_prefix("$(prefix)_$(freq_off)", 0:9))
end

const all_imgs = [load_images(joinpath(@__DIR__, "data", "BaOven_20210730_current"),
                              currents[i], times[i]) for i in 1:length(currents)]

const prefix = joinpath(@__DIR__, "imgs", "BaOven_current_20210730")

const mean_imgs = [(mean(imgs.res) .- mean(imgs.off)) ./ time
                   for (time, imgs) in zip(times, all_imgs)]

const signal_roi = (115:131, 135:162)

const ncurrents = length(currents)
const nrows = ceil(Int, sqrt(ncurrents))
const ncols = ceil(Int, ncurrents / nrows)

const mean_imgs_max = maximum(maximum(img[signal_roi...]) for img in mean_imgs)
const mean_imgs_min = minimum(minimum(img[signal_roi...]) for img in mean_imgs)

figure(figsize=[5.4, 4.8] .* 1.2)
for i in 1:ncurrents
    subplot(nrows, ncols, i)
    imshow(mean_imgs[i][signal_roi...], vmin=mean_imgs_min, vmax=mean_imgs_max)
    title("$(currents[i]) A", fontsize=12)
    xticks([])
    yticks([])
end
tight_layout(pad=0.3)
NaCsPlot.maybe_save("$(prefix)_imgs_norm")

figure(figsize=[5.4, 4.8] .* 1.2)
for i in 1:ncurrents
    subplot(nrows, ncols, i)
    imshow(mean_imgs[i][signal_roi...])
    title("$(currents[i]) A", fontsize=12)
    xticks([])
    yticks([])
end
tight_layout(pad=0.3)
NaCsPlot.maybe_save("$(prefix)_imgs_auto")

function get_single_count_per_sec(img, time)
    return sum(@view img[signal_roi...]) / time
end

function get_count_per_sec(imgs, time)
    res = Float64[]
    off = Float64[]
    for img in imgs.res
        push!(res, get_single_count_per_sec(img, time))
    end
    for img in imgs.off
        push!(off, get_single_count_per_sec(img, time))
    end
    return Unc(mean(res) - mean(off),
               sqrt(std(res)^2 / length(res) + std(off)^2 / length(off)))
end

const count_per_sec_uncs = [get_count_per_sec(all_imgs[i], times[i])
                            for i in 1:length(currents)]
# @show count_per_sec_uncs

function uncs_to_errorbar(uncs)
    as = Float64[]
    ss = Float64[]
    for unc in uncs
        push!(as, unc.a)
        push!(ss, unc.s)
    end
    return as, ss
end

# Objective focal length: 100mm, clear diameter: 0.9inch (from retaining ring) / 22.86mm
# Imaging lens focal length: 175mm
# Pixel size 6.5um (100% filling AFAICT)
# Quantum efficiency: 78%
# Factor from read mode: 16 ??
# Binning: 8
# Beam size: e^-2 FW 160um
# Viewing angle: 45deg
# Atom beam width: ~580um
# Florescence area size: 92800um^2 -> 0.000928cm^2

const QE = 0.78 / 16
const transmission = 1
const area = 0.000928 # cm^2
const solid_angle_ratio = (100 - sqrt(100^2 - (22.86 / 2)^2)) / (2 * 100) # h / 2R
const m_Ba = 138e-3 / 6.02e23
const k_B = 1.38e-23
const T_Ba = 1000
const v_Ba = sqrt(k_B * T_Ba / m_Ba) * 100 # in cm/s

const atom_per_cm2_per_sec_uncs = count_per_sec_uncs ./ QE ./ solid_angle_ratio ./ area ./ transmission
const atom_per_cm3_uncs = atom_per_cm2_per_sec_uncs ./ v_Ba

figure()
errorbar(currents, uncs_to_errorbar(atom_per_cm2_per_sec_uncs)...)
grid()
yscale("log")
xlabel("Oven current (A)")
ylabel("Flux (cm\$^{-2}\\cdot s^{-1}\$)")
title("Ba Flux")
NaCsPlot.maybe_save("$(prefix)_flux")

figure()
errorbar(currents, uncs_to_errorbar(atom_per_cm3_uncs)...)
grid()
yscale("log")
xlabel("Oven current (A)")
ylabel("Density (cm\$^{-3}\$)")
title("Ba Density (assume $(T_Ba) K)")
NaCsPlot.maybe_save("$(prefix)_density")

NaCsPlot.maybe_show()

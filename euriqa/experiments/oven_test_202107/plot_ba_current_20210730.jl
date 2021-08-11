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

NaCsPlot.maybe_show()

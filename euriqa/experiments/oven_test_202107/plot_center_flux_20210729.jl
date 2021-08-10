#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
import NaCsCalc.Format: Unc, Sci
using Images
using Statistics

const off_names = ["3s_0A_261uW_1.sif",
                   "3s_0A_261uW_2.sif",
                   "3s_0A_261uW_3.sif",
                   "3s_0A_261uW_4.sif",
                   "3s_0A_261uW_5.sif",
                   "3s_0A_261uW_6.sif",
                   "3s_0A_261uW_7.sif",
                   "3s_0A_261uW_8.sif",
                   "3s_0A_261uW_9.sif"]
const on_names = ["3s_2p25A_261uW_1.sif",
                  "3s_2p25A_261uW_2.sif",
                  "3s_2p25A_261uW_3.sif",
                  "3s_2p25A_261uW_4.sif",
                  "3s_2p25A_261uW_5.sif",
                  "3s_2p25A_261uW_6.sif",
                  "3s_2p25A_261uW_7.sif",
                  "3s_2p25A_261uW_8.sif",
                  "3s_2p25A_261uW_9.sif"]

const off_imgs = [Float64.(load(joinpath(@__DIR__, "data",
                                         "YbOven_20210729_center", name)).data[65:143, 81:159, 1])
                  for name in off_names]
const on_imgs = [Float64.(load(joinpath(@__DIR__, "data",
                                        "YbOven_20210729_center", name)).data[65:143, 81:159, 1])
                 for name in on_names]

const prefix = joinpath(@__DIR__, "imgs", "YbOven_center_20210729")

const off_img_avg = mean(off_imgs)
const on_img_avg = mean(on_imgs)

figure()
imshow(off_img_avg, vmax=300)
colorbar()
xlabel("x (pixel)")
ylabel("y (pixel)")
title("Oven off")
NaCsPlot.maybe_save("$(prefix)_off")

figure()
imshow(on_img_avg, vmax=300)
colorbar()
xlabel("x (pixel)")
ylabel("y (pixel)")
title("Oven on")
NaCsPlot.maybe_save("$(prefix)_on")

figure()
imshow(on_img_avg .- off_img_avg, vmin=0)
colorbar()
xlabel("x (pixel)")
ylabel("y (pixel)")
title("Difference")
NaCsPlot.maybe_save("$(prefix)_diff")

NaCsPlot.maybe_show()

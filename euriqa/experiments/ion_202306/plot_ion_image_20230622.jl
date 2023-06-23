#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Format: Unc
using NaCsData
using NaCsPlot
using PyPlot
using MAT
using Statistics

read_mat_image(fname) = matopen(fname) do mat
    data = read(mat, "data")
    return data[:, :]
end

function average_image(images)
    a = zeros(size(images[1]))
    s = zeros(size(images[1]))
    for img in images
        a .= a .+ img
        s .= s .+ img.^2
    end
    nimgs = length(images)
    a .= a ./ nimgs
    s .= s ./ nimgs
    s .= (s .- a.^2) ./ (nimgs - 1)
    return a, s
end

const ion_img, ion_img_s2 = average_image([read_mat_image(joinpath(@__DIR__, "data",
                                                                   "1ion_$i.mat"))
                                           for i in 1:10])
const bg_img, bg_img_s2 = average_image([read_mat_image(joinpath(@__DIR__, "data",
                                                                 "background_$i.mat"))
                                         for i in 1:10])

const diff_img = ion_img .- bg_img
const diff_img_s2 = ion_img_s2 .+ bg_img_s2

function integrate_circle(img, x0, y0, r)
    # using a very simple edge detection/anti-aliasing algorithm
    x1 = max(floor(Int, x0 - r), 1)
    x2 = min(ceil(Int, x0 + r), size(img, 1))
    y1 = max(floor(Int, y0 - r), 1)
    y2 = min(ceil(Int, y0 + r), size(img, 2))

    s = 0.0
    for x in x1:x2
        for y in y1:y2
            v = img[x, y]
            d = hypot(x - x0, y - y0)
            if d <= r - 0.5
                s += v
            elseif d <= r + 0.5
                s += v * (r + 0.5 - d)
            end
        end
    end
    return s
end

const prefix = joinpath(@__DIR__, "imgs", "ion_image_20230622")

const r_fiber = 9
const d_fiber = 19

const center = integrate_circle(diff_img, 784, 1520, r_fiber)
const center_s2 = integrate_circle(diff_img_s2, 784, 1520, r_fiber)

const angles = range(0, 2π, 181)

function compute_crosstalk(angles, r_fiber, d_fiber)
    side_count = [integrate_circle(diff_img,
                                   784 + d_fiber * cos(angle),
                                   1520 + d_fiber * sin(angle), r_fiber)
                  for angle in angles]
    side_count_s2 = [integrate_circle(diff_img_s2,
                                      784 + d_fiber * cos(angle),
                                      1520 + d_fiber * sin(angle), r_fiber)
                     for angle in angles]
    xc = side_count ./ center
    xc_s = sqrt.(side_count_s2 ./ side_count.^2 .+ center_s2 / center.^2) .* abs.(xc)
    return xc, xc_s
end

const xc1, xc1_s = compute_crosstalk(angles, r_fiber, d_fiber)
const xc2, xc2_s = compute_crosstalk(angles, r_fiber, d_fiber * 2)
const xc3, xc3_s = compute_crosstalk(angles, r_fiber, d_fiber * 3)
const xc5, xc5_s = compute_crosstalk(angles, r_fiber, d_fiber * 5)
const xc8, xc8_s = compute_crosstalk(angles, r_fiber, d_fiber * 8)
const xc13, xc13_s = compute_crosstalk(angles, r_fiber, d_fiber * 13)

figure()
errorbar(angles .* (180 / π), xc1 .* 100, xc1_s .* 100, label="dist 1")
errorbar(angles .* (180 / π), xc2 .* 100, xc2_s .* 100, label="dist 2")
errorbar(angles .* (180 / π), xc3 .* 100, xc3_s .* 100, label="dist 3")
errorbar(angles .* (180 / π), xc5 .* 100, xc5_s .* 100, label="dist 5")
errorbar(angles .* (180 / π), xc8 .* 100, xc8_s .* 100, label="dist 8")
errorbar(angles .* (180 / π), xc13 .* 100, xc13_s .* 100, label="dist 13")
legend(fontsize=13, ncol=3)
grid()
xlabel("Direction (\$^\\circ\$)")
ylabel("Crosstalk (%)")
xlim([0, 360])
ylim([0, 8])
NaCsPlot.maybe_save("$(prefix)_crosstalk")

NaCsPlot.maybe_show()

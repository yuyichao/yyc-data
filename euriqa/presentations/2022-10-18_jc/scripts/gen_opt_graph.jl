#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using PyPlot
using NaCsCalc
using NaCsPlot

function gen_center_gaussian(size, x2, xy, y2)
    r = range(-1, 1, size)
    return [exp(-(x2 * x^2 / 2 + xy * x * y + y2 * y^2 / 2)) for x in r, y in r]
end

const image_size = 1024
const sym_image = gen_center_gaussian(image_size, 5, 0, 5)
const asym_image = gen_center_gaussian(image_size, 18, 17, 18)

val_to_idx(x) = (x + 1) * image_size / 2 - 0.5

const prefix = joinpath(@__DIR__, "../imgs/opt")

figure()
gca().axis("off")
imshow(sym_image, cmap="YlGn")
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_sym", pad_inches=0)

figure()
gca().axis("off")
imshow(sym_image, cmap="YlGn")
plot(val_to_idx(0.8), val_to_idx(-0.7), "o", color="red")
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_sym_point", pad_inches=0)

figure()
gca().axis("off")
imshow(sym_image, cmap="YlGn")
plot([val_to_idx(0.8), val_to_idx(0)], [val_to_idx(-0.7), val_to_idx(0)],
     color="salmon")
plot([val_to_idx(0.8), val_to_idx(0), val_to_idx(0)],
     [val_to_idx(-0.7), val_to_idx(-0.7), val_to_idx(0)],
     color="orange")
plot(val_to_idx(0.8), val_to_idx(-0.7), "o", color="red")
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_sym_align", pad_inches=0)

figure()
gca().axis("off")
imshow(asym_image, cmap="YlGn")
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_asym", pad_inches=0)

figure()
gca().axis("off")
imshow(asym_image, cmap="YlGn")
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(-1), val_to_idx(1)],
     "--", color="gray", alpha=0.5)
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(1), val_to_idx(-1)],
     "--", color="gray", alpha=0.5)
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_asym_axis", pad_inches=0)

figure()
gca().axis("off")
imshow(asym_image, cmap="YlGn")
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(-1), val_to_idx(1)],
     "--", color="gray", alpha=0.5)
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(1), val_to_idx(-1)],
     "--", color="gray", alpha=0.5)
plot(val_to_idx(0.8), val_to_idx(-0.7), "o", color="red")
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_asym_point", pad_inches=0)

const naive_path = [(0.8, -0.7)]

for i in 1:30
    p = naive_path[end]
    if abs(p[1]) < abs(p[2])
        push!(naive_path, (p[1], -p[1] * 8.5/9))
    else
        push!(naive_path, (-p[2] * 8.5/9, p[2]))
    end
end

figure()
gca().axis("off")
imshow(asym_image, cmap="YlGn")
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(-1), val_to_idx(1)],
     "--", color="gray", alpha=0.5)
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(1), val_to_idx(-1)],
     "--", color="gray", alpha=0.5)
plot([val_to_idx(p[1]) for p in naive_path], [val_to_idx(p[2]) for p in naive_path],
     color="orange")
plot(val_to_idx(0.8), val_to_idx(-0.7), "o", color="red")
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_asym_naive", pad_inches=0)

const jump_path = [(0.8, -0.7)]

for i in 1:10
    p = jump_path[end]
    if abs(p[1]) < abs(p[2])
        new_p2 = -p[1] * 5/9
        push!(jump_path, (p[1], new_p2))
    else
        new_p1 = -p[2] * 5/9
        push!(jump_path, (new_p1, p[2]))
    end
end

figure()
gca().axis("off")
imshow(asym_image, cmap="YlGn")
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(-1), val_to_idx(1)],
     "--", color="gray", alpha=0.5)
plot([val_to_idx(1), val_to_idx(-1)], [val_to_idx(1), val_to_idx(-1)],
     "--", color="gray", alpha=0.5)
plot([val_to_idx(p[1]) for p in naive_path], [val_to_idx(p[2]) for p in naive_path],
     color="orange")
plot([val_to_idx(p[1]) for p in jump_path], [val_to_idx(p[2]) for p in jump_path],
     color="salmon")
plot(val_to_idx(0.8), val_to_idx(-0.7), "o", color="red")
xticks([])
yticks([])
tick_params(which="both", length=0, width=0, pad=0)
NaCsPlot.maybe_save("$(prefix)_asym_jump", pad_inches=0)

NaCsPlot.maybe_show()

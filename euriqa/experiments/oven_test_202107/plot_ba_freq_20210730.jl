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

function load_sif_prefix(prefix, idxs)
    return [Float64.(load("$(prefix)_$(idx).sif").data) for idx in idxs]
end

const freqs = [725258630, 725258650, 725258670, 725258690, 725258710, 725258730,
               725258750, 725258770, 725258790, 725258800, 725258810, 725258820,
               725258829, 725258830, 725258840, 725258850, 725258860, 725258870,
               725258880, 725258900, 725258920, 725258940, 725258960, 725258980,
               725259000, 725259020, 725259040, 725259060, 725259080, 725259100, 725259120]

const all_imgs = [load_sif_prefix(joinpath(@__DIR__, "data", "BaOven_20210730_freq",
                                       "3s_3p8A_754uW_$freq"), 0:9)[:, :, 1]
              for freq in freqs]

const prefix = joinpath(@__DIR__, "imgs", "BaOven_freq_20210730")

const mean_imgs = [mean(imgs) for imgs in all_imgs]

const signal_roi = (115:131, 135:165)
const background_roi = (90:156, 100:200)

const nfreqs = length(freqs)
const nrows = ceil(Int, sqrt(nfreqs))

const mean_imgs_max = maximum(maximum(img[signal_roi...]) for img in mean_imgs)
const mean_imgs_min = minimum(minimum(img[signal_roi...]) for img in mean_imgs)

figure(figsize=[6.4, 4.8] .* 2)
for i in 1:nfreqs
    subplot(nrows, nrows, i)
    imshow(mean_imgs[i][signal_roi...], vmin=mean_imgs_min, vmax=mean_imgs_max)
    title("$(freqs[i]) MHz", fontsize=12)
    xticks([])
    yticks([])
end
tight_layout(pad=0.3)
NaCsPlot.maybe_save("$(prefix)_imgs_norm")

figure(figsize=[6.4, 4.8] .* 2)
for i in 1:nfreqs
    subplot(nrows, nrows, i)
    imshow(mean_imgs[i][signal_roi...])
    title("$(freqs[i]) MHz", fontsize=12)
    xticks([])
    yticks([])
end
tight_layout(pad=0.3)
NaCsPlot.maybe_save("$(prefix)_imgs_auto")

function get_single_signal(img)
    bg_img = @view img[background_roi...]
    sig_img = @view img[signal_roi...]

    bg_mean = (sum(bg_img) - sum(sig_img)) / (length(bg_img) - length(sig_img))
    return sum(sig_img) .- bg_mean * length(sig_img)
end

const signals = Float64[]
const signal_uncs = Float64[]

for imgs in all_imgs
    ss = [get_single_signal(img) for img in imgs]
    push!(signals, mean(ss))
    push!(signal_uncs, std(ss) / sqrt(length(ss)))
end

function model(x, p)
    return p[1] .+ p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2) .+ (x .- p[3]) .* p[5]
end
function model2(x, p)
    return p[1] .+ p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2) .+ p[5] ./ (1 .+ ((x .- p[6]) ./ (p[7] / 2)).^2)
end

fit = fit_data(model, freqs[1:24], signals[1:24],
               [0.0, 8000, 725258829, 50, 0.0])
fit2 = fit_data(model2, freqs, signals, signal_uncs, [fit.param[1:4]; 1300; 725259060; 100])
@show fit2.param

figure()
errorbar(freqs .- fit2.param[3], signals, signal_uncs, fmt=".", color="C0")
plot(fit2.plotx .- fit2.param[3], fit2.ploty, color="C0")
text(20, 8000, "\$f_0=$(round(Int, fit2.param[3]))\$ MHz", fontsize=15, color="C0")
text(60, 7300, "\$\\Gamma_0=$(round(Int, fit2.param[4]))\$ MHz", fontsize=15, color="C0")
text(63, 2000, "\$f'=$(round(Int, fit2.param[6]))\$ MHz", fontsize=15, color="C0")
text(100, 2700, "\$\\Gamma'=$(round(Int, fit2.param[7]))\$ MHz", fontsize=15, color="C0")
grid()
ylim([0, ylim()[2]])
xlabel("Detuning (MHz)")
ylabel("Total Count")
NaCsPlot.maybe_save("$(prefix)_fit")

NaCsPlot.maybe_show()

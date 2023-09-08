#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot

function load_data(name)
    data = readdlm(joinpath(@__DIR__, "data", "$(name).csv"),
                   ',', skipstart=2)
    for i in 1:size(data, 1)
        if any(x->!isa(x, Float64), @view(data[i, :]))
            return data[1:(i - 1), :]
        end
    end
    return data
end

const prefix = joinpath(@__DIR__, "imgs/scope_trace_20230907")

const data1 = load_data("871_Cavity_Signal_20230907_1")
const data2 = load_data("871_Cavity_Signal_20230907_2")

figure()
plot(data1[:, 1] .* 1e6, data1[:, 2] ./ maximum(abs.(data1[:, 2])),
     label="Error signal")
plot(data1[:, 1] .* 1e6, data1[:, 4] ./ maximum(abs.(data1[:, 4])),
     label="Transmission")
grid()
xlabel("Time (\$\\mu s\$)")
ylabel("(arb. unit)")
legend(fontsize=14)
NaCsPlot.maybe_save("$(prefix)1")

figure()
plot(data2[:, 1] .* 1e6, data2[:, 2] ./ maximum(abs.(data2[:, 2])),
     label="Error signal")
plot(data2[:, 1] .* 1e6, data2[:, 4] ./ maximum(abs.(data2[:, 4])),
     label="Transmission")
grid()
xlabel("Time (\$\\mu s\$)")
ylabel("(arb. unit)")
legend(fontsize=14)
NaCsPlot.maybe_save("$(prefix)2")

NaCsPlot.maybe_show()

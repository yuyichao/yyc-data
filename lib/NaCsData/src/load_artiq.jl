#!/usr/bin/julia -f

using HDF5

function load_dax_scan1(fname::AbstractString)
    h5open(load_dax_scan1, fname)
end

# 1D scan, single parameter, single measurement, single ion
function load_dax_scan1(fd)
    thresh = read(fd, "archive/system.pmt.state_detection_threshold")
    scan = delete!(read(fd, "datasets/scan"), "product")
    @assert length(scan) == 1
    param_name = first(keys(scan))
    params = first(values(scan))
    nparams = length(params)
    pmt_counts = read(fd, "datasets/histogram_context/histogram/raw")
    actual_nparams = length(pmt_counts)
    if nparams == actual_nparams
        counts = Matrix{Int}(undef, nparams, 2)
        for i in 1:nparams
            pmt_count = pmt_counts[string(i - 1)]
            counts[i, 1] = length(pmt_count)
            counts[i, 2] = count(>=(thresh), pmt_count)
        end
    else
        counts = Matrix{Int}(undef, actual_nparams, 2)
        new_params = similar(params, nparams)
        idx = 0
        for i in 1:nparams
            k = string(i - 1)
            (k in keys(pmt_counts)) || continue
            idx += 1
            new_params[idx] = params[i]
            pmt_count = pmt_counts[k]
            counts[idx, 1] = length(pmt_count)
            counts[idx, 2] = count(>=(thresh), pmt_count)
        end
        params = new_params
    end
    return (param_name,), CountData(params, counts)
end

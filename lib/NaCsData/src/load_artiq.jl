#!/usr/bin/julia -f

using HDF5

function load_dax_scan1(fname::AbstractString; thresh=nothing)
    h5open(fh->load_dax_scan1(fh, thresh=thresh), fname)
end

function _convert_count(pmt_counts, thresh, si=1)
    nparams = length(pmt_counts)
    counts = Matrix{Int}(undef, nparams, 2)
    for i in 1:nparams
        pmt_count = @view pmt_counts[i][si, :]
        counts[i, 1] = length(pmt_count)
        counts[i, 2] = count(>=(thresh), pmt_count)
    end
    return counts
end

# 1D scan, single parameter, single measurement, single ion
function load_dax_scan1(fd; thresh=nothing)
    if thresh === nothing
        thresh = read(fd, "archive/system.pmt.state_detection_threshold")
    end
    scan = delete!(read(fd, "datasets/scan"), "product")
    @assert length(scan) == 1
    param_name = first(keys(scan))
    params = first(values(scan))
    nparams = length(params)
    pmt_counts = read(fd, "datasets/histogram_context/histogram/raw")
    actual_nparams = length(pmt_counts)
    pmt_counts_ary = Vector{Matrix{Int}}(undef, actual_nparams)
    if nparams == actual_nparams
        for i in 1:nparams
            pmt_counts_ary[i] = pmt_counts[string(i - 1)]
        end
    else
        new_params = similar(params, actual_nparams)
        idx = 0
        for i in 1:nparams
            k = string(i - 1)
            (k in keys(pmt_counts)) || continue
            idx += 1
            new_params[idx] = params[i]
            pmt_counts_ary[idx] = pmt_counts[k]
        end
        params = new_params
        nparams = actual_nparams
    end
    @assert nparams > 0
    nsites = size(pmt_counts_ary[1], 1)
    return (param_name,), [CountData(params, _convert_count(pmt_counts_ary,
                                                            thresh, si))
                           for si in 1:nsites]
end

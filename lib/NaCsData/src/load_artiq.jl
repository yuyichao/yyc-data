#!/usr/bin/julia -f

using HDF5

function load_dax_scan_logicals1(fname::AbstractString; kws...)
    h5open(fh->load_dax_scan_logicals1(fh; kws...), fname)
end

# Return nshots * nsites * nseqs
function _conert_count_to_logical(params, pmt_counts, thresh)
    nparams = length(pmt_counts)
    nseqs = 0
    nsites = 0
    for i in 1:nparams
        pmt_count = pmt_counts[i]
        _nsites = size(pmt_count, 1)
        @assert _nsites != 0
        if nsites == 0
            nsites = _nsites
        else
            @assert nsites == _nsites
        end
        nseqs += size(pmt_count, 2)
    end
    logicals = Array{Bool,3}(undef, 1, nsites, nseqs)
    seq_params = Vector{eltype(params)}(undef, nseqs)
    seq_offset = 0
    for i in 1:nparams
        pmt_count = pmt_counts[i]
        nrep = size(pmt_count, 2)
        seq_params[seq_offset .+ (1:nrep)] .= Ref(params[i])
        logicals[1, :, seq_offset .+ (1:nrep)] .= pmt_count[:, 1:nrep] .> thresh
        seq_offset += nrep
    end
    return seq_params, logicals
end

# 1D scan, single parameter, single measurement
function load_dax_scan_logicals1(fd; thresh=nothing, index_param=false)
    thresh, param_names, params, pmt_counts_ary, nsites =
        _load_dax_scan(Int, fd, thresh, index_param, nothing)
    return param_names, _conert_count_to_logical(params, pmt_counts_ary, thresh)
end

function load_dax_scan1(fname::AbstractString; kws...)
    h5open(fh->load_dax_scan1(fh; kws...), fname)
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

function _load_dax_scan(::Type{T}, fd, thresh, index_param, count_processor) where T
    if thresh === nothing
        thresh = read(fd, "archive/system.pmt.state_detection_threshold")
    end
    scan = read(fd, "datasets/scan/product")
    param_names = (sort(keys(scan))...,)
    param_values = getindex.(Ref(scan), param_names)
    nparams = length(param_values[1])
    if index_param
        param_names = ("Index",)
        params = collect(1:nparams)
    elseif length(param_names) == 1
        # Backward compatible
        params = param_values[1]
    else
        names = Symbol.(param_names)
        NT = NamedTuple{names}
        params = [NT(getindex.(param_values, i)) for i in 1:nparams]
    end
    pmt_counts = read(fd, "datasets/histogram_context/histogram/raw")
    actual_nparams = length(pmt_counts)
    pmt_counts_ary = Vector{Matrix{T}}(undef, actual_nparams)
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
    if count_processor !== nothing
        for pmt_counts in pmt_counts_ary
            count_processor(pmt_counts)
        end
    end
    nsites = size(pmt_counts_ary[1], 1)

    return thresh, param_names, params, pmt_counts_ary, nsites
end

# 1D scan, single parameter, single measurement
function load_dax_scan1(fd; thresh=nothing, count_processor=nothing, index_param=false)
    thresh, param_names, params, pmt_counts_ary, nsites =
        _load_dax_scan(Float64, fd, thresh, index_param, count_processor)
    return param_names, [CountData(params, _convert_count(pmt_counts_ary,
                                                          thresh, si))
                         for si in 1:nsites]
end

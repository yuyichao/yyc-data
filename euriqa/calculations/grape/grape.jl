#!/usr/bin/julia

abstract type AbstractOP{VT,NParams} end

function compute end
function set_params end

(nparams(::Type{T} where T<:AbstractOP{VT,NParams}) where {VT,NParams}) = NParams

_nparams_ops(params_types) = sum(nparams, params_types)

struct Sequence{VT,OPS<:NTuple{N,AbstractOP} where N}
    ops::OPS
    val_buff::Vector{VT}
    grad_buff::Vector{VT}
    prefix_buff::Vector{VT}
    suffix_buff::Vector{VT}

    function Sequence{VT}(ops::OPS) where OPS<:NTuple{N,AbstractOP{VT}} where {VT,N}
        @assert N > 0
        return new{VT,OPS}(ops, Vector{VT}(undef, N),
                           Vector{VT}(undef, _nparams_ops(OPS.parameters)),
                           Vector{VT}(undef, N - 1), Vector{VT}(undef, N - 1))
    end
end

@generated function nparams(::Type{Sequence{VT,OPS}}) where {VT,OPS}
    return _nparams_ops(OPS.parameters)
end

@generated function param_range(::Type{Sequence{VT,OPS}}) where {VT,OPS}
    OP_Types = (OPS.parameters...,)
    op_nparams = nparams.(OP_Types)
    nops = length(OP_Types)
    cum_nparams = cumsum(op_nparams)
    starts = (0, cum_nparams[1:end - 1]...) .+ 1
    return starts, cum_nparams
end

@generated function set_params(s::Sequence, params)
    ex = quote
        @assert length(params) == nparams($s)
        ops = s.ops
    end
    starts, ends = param_range(s)
    for (i, (start_idx, end_idx)) in enumerate(zip(starts, ends))
        push!(ex.args, :(@inline set_params(@inbounds(ops[$i]),
                                            @view params[$start_idx:$end_idx])))
    end
    push!(ex.args, :(return))
    return ex
end

@generated function _eval_compute(s)
    ex = quote
        ops = s.ops
        grad_buff = s.grad_buff
        val_buff = s.val_buff
    end
    starts, ends = param_range(s)
    for (i, (start_idx, end_idx)) in enumerate(zip(starts, ends))
        push!(ex.args, :(
            @inbounds begin
                @inline val_buff[$i] =
                    compute(ops[$i], @view grad_buff[$start_idx:$end_idx])
            end))
    end
    push!(ex.args, :(return))
    return ex
end

function compute(s::Sequence, grads)
    @assert length(grads) == nparams(typeof(s))
    N = length(s.ops)
    starts, ends = param_range(typeof(s))
    @assert N > 0
    if N == 1
        return compute(s.ops[1], grads)
    end
    @inline _eval_compute(s)
    @inbounds s.prefix_buff[1] = s.val_buff[1]
    @inbounds for i in 2:N - 1
        s.prefix_buff[i] = s.prefix_buff[i - 1] * s.val_buff[i]
    end
    res = @inbounds s.prefix_buff[N - 1] * s.val_buff[N]

    @inbounds s.suffix_buff[N - 1] = s.val_buff[N]
    @inbounds for i in N - 2:-1:1
        s.suffix_buff[i] = s.val_buff[i + 1] * s.suffix_buff[i + 1]
    end
    @inbounds for op_idx in 1:N
        pstart = starts[op_idx]
        pend = ends[op_idx]
        if op_idx == 1
            suffix = s.suffix_buff[op_idx]
            for param_idx in pstart:pend
                grads[param_idx] = s.grad_buff[param_idx] * suffix
            end
        elseif op_idx == N
            prefix = s.prefix_buff[op_idx - 1]
            for param_idx in pstart:pend
                grads[param_idx] = prefix * s.grad_buff[param_idx]
            end
        else
            prefix = s.prefix_buff[op_idx - 1]
            suffix = s.suffix_buff[op_idx]
            for param_idx in pstart:pend
                grads[param_idx] = prefix * s.grad_buff[param_idx] * suffix
            end
        end
    end
    return res
end

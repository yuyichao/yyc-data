#!/usr/bin/julia

function average_data(datas; offset=0)
    npoints = sum(length(data) for data in datas)
    total = sum(sum(data) for data in datas)
    return total / npoints + offset
end

function integrate!(out, data; offset=0)
    T = eltype(out)
    v0 = zero(T)
    out[1] = v0
    d0 = data[1]
    for i in 2:length(data)
        d1 = data[i]
        v1 = v0 + T(d0 + d1) / 2 + offset
        out[i] = v1
        d0 = d1
        v0 = v1
    end
    return out
end

function get_integrated(int_data, data, idx; offset=0)
    T = eltype(int_data)
    I = floor(Int, idx)
    v0 = int_data[I]
    if I == idx
        return v0
    end
    d0 = data[I]
    d1 = data[I + 1]
    α = idx - I
    return T(evalpoly(idx - I, (v0, d0 + offset, (d1 - d0) / 2)))
end

function calc_phase!(out, data, flip; offset=0)
    integrate!(out, data; offset=offset)
    if flip == 0
        return out
    end
    get_t(t) = get_integrated(out, data, t + 1; offset=offset)
    @inbounds for t in length(out) - 1:-1:0
        v = out[t + 1]
        sign = 1
        for j in flip:-1:1
            sign = -sign
            v += sign * 2 * get_t(j * t / (flip + 1))
        end
        out[t + 1] = v
    end
    return out
end

function calc_rabi!(out, data, flip, buf, weight, Ω; offset=0)
    calc_phase!(buf, data, flip; offset=offset)
    out .= muladd.(muladd.(-1, cos.(buf .* Ω), 1), weight, out)
    return out
end

function sample_rabi(datas, len, num, flip, Ω; offset=0, yscale=1)
    T = eltype(datas[1])
    buf = Vector{T}(undef, len)
    out = zeros(T, len)
    for _ in 1:num
        data = rand(datas)
        dlen = length(data)
        @assert dlen > len * 2
        starti = rand(1:(dlen - len + 1))
        calc_rabi!(out, @view(data[starti:starti + len - 1]), flip,
                   buf, 0.5 * yscale / num, Ω; offset=offset)
    end
    return out
end

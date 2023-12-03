#!/usr/bin/julia

struct Params{T,D}
    P_turn::T # 1 - exp(με) or 1 - exp(-με)
    P_hop::T # tε
    forward_turn::Bool # 1 > exp(με)
    ntime::Int
    nspace::NTuple{D,Int}
    function Params{T}(μ, ε, t, ntime, nspace::NTuple{D}) where {T,D}
        με = μ * ε
        return new{T,D}(1 - exp(-abs(με)), t * ε, με < 0, ntime, nspace)
    end
end

struct State{T,D,_D}
    param::Params{T,D}
    graph::Array{UInt8,_D}
    function State(param::Params{T,D}) where {T,D}
        _D = D + 1
        return new{T,D,_D}(param, zeros(UInt8, (param.ntime, param.nspace...)))
    end
end

# 0: empty
# 1: not moving
# 2d: d'th dimension, -1
# 2d + 1: d'th dimension, +1
# First element is the direction forward, second element is direction backward
@inline function get_dir(state::State, idx)
    v = @inbounds state.graph[idx...]
    return (v & 0xf), (v >> 4)
end

@inline function set_dir(state::State, idx, dirs)
    @inbounds state.graph[idx...] = (dirs[1] | (dirs[2] << 4)) % UInt8
    return
end

@inline function get_next(state::State{T,D}, pidx, dir) where {T,D}
    # @assert dir != 0
    if dir == 1
        return pidx
    end
    _dir = dir >> 1
    # @assert 1 <= _dir <= D
    nspace = state.param.nspace
    if dir & 1 == 0
        _pidx = ntuple(Val(D)) do i
            v = @inbounds pidx[i]
            if i != _dir
                return v
            end
            v += 1
            @inbounds if v > nspace[i]
                return 1
            end
            return v
        end
    else
        _pidx = ntuple(Val(D)) do i
            v = @inbounds pidx[i]
            if i != _dir
                return v
            end
            v -= 1
            @inbounds if v <= 0
                return nspace[i]
            end
            return v
        end
    end
    # @assert pidx != _pidx
    return _pidx
end

function gather_result(state)
    nstraight = 0
    nhop = 0
    nlines = 0
    @inbounds for si in CartesianIndices(state.param.nspace)
        dir0, = get_dir(state, (1, si.I...))
        if dir0 == 0
            continue
        end
        nlines += 1
        if dir0 == 1
            nstraight += 1
        else
            nhop += 1
        end
        pidx = get_next(state, si.I, dir0)
        for ti in 2:state.param.ntime
            dir, = get_dir(state, (ti, pidx...))
            if dir == 1
                nstraight += 1
            else
                nhop += 1
            end
            pidx = get_next(state, pidx, dir)
        end
    end
    return nlines, nstraight, nhop
end

struct Step{D}
    forward::Bool
    tidx::Int
    pidx::NTuple{D,Int}

    # For forward direction moving, this is the backward direction to be written
    # to the new site.
    # For backward direction moving, this is the direction to move
    dirb::UInt8
end

# For forward step, this function will write to the current site
# For backward step, this function will assume that the current site
# is written to already
@inline function sweep_iterate!(state::State{T,D}, step::Step{D}) where {T,D}
    param = state.param
    if step.forward
        # Next time
        tidx_next = step.tidx == param.ntime ? 1 : step.tidx + 1

        # Next coordinate
        r = floor(Int, rand(T) / param.P_hop)
        if r >= 2 * D
            dirf = 0x1
            dirb_next_new = 0x1
        else
            dirf = (r + 2) % UInt8
            dirb_next_new = dirf ⊻ 0x1
        end
        set_dir(state, (step.tidx, step.pidx...), (dirf, step.dirb))
        pidx_next = get_next(state, step.pidx, dirf)
        dirf_next, dirb_next_orig = get_dir(state, (tidx_next, pidx_next...))
        if dirf_next == 0
            # Empty site
            return false, Step{D}(!param.forward_turn || rand(T) >= param.P_turn,
                                  tidx_next, pidx_next, dirb_next_new)
        end
        set_dir(state, (tidx_next, pidx_next...), (dirf_next, dirb_next_new))
        if dirb_next_orig == 0
            # Found a start, finish!
            return true, Step{D}(false, tidx_next, pidx_next, 0x0)
        end
        # Found an intersection
        return false, Step{D}(false, tidx_next, pidx_next, dirb_next_orig)
    else
        # Next time
        tidx_next = step.tidx - 1
        if tidx_next <= 0
            tidx_next = param.ntime
        end

        # Next coordinate
        pidx_next = get_next(state, step.pidx, step.dirb)
        dirf_next, dirb_next = get_dir(state, (tidx_next, pidx_next...))
        set_dir(state, (tidx_next, pidx_next...), (0x0, 0x0))

        if !param.forward_turn && rand(T) < param.P_turn
            # Turn around
            return false, Step{D}(true, tidx_next, pidx_next, dirb_next)
        end
        return dirb_next == 0, Step{D}(false, tidx_next, pidx_next, dirb_next)
    end
end

function sweep!(state::State{T,D}) where {T,D}
    param = state.param
    tidx, pidx... = rand(CartesianIndices((param.ntime, param.nspace...))).I
    dirf, dirb = get_dir(state, (tidx, pidx...))

    step = Step{D}(dirb == 0, tidx, pidx, dirb)

    while true
        finish, step = sweep_iterate!(state, step)
        if finish
            break
        end
    end

    return gather_result(state)
end

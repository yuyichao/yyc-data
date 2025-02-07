#!/usr/bin/julia

struct Params{T,D}
    P_turn::T # 1 - exp(με) or 1 - exp(-με)
    P_hop::T # tε
    forward_turn::Bool # 1 > exp(με)
    ntime::Int
    nspace::NTuple{D,Int}
    Estraight::T
    Ehop::T
    function Params{T}(μ, ε, t, ntime, nspace::NTuple{D}) where {T,D}
        με = μ * ε
        return new{T,D}(1 - exp(-abs(με)), t * ε, με < 0, ntime, nspace,
                        2 * D * t / (1 - 2 * D * t * ε) / ntime, -1 / ε / ntime)
    end
end

mutable struct State{T,D,_D}
    const param::Params{T,D}
    const graph::Array{UInt8,_D}
    nstraight::Int
    nhop::Int
    nlines::Int
    winding::Int
    χ₀::Int
    function State(param::Params{T,D}) where {T,D}
        @assert D <= 7 # The UInt8 state can support up to 7 dimensions
        _D = D + 1
        return new{T,D,_D}(param, zeros(UInt8, (param.ntime, param.nspace...)),
                           0, 0, 0, 0, 0)
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
        return pidx, 0
    end
    _dir = dir >> 1
    # @assert 1 <= _dir <= D
    nspace = state.param.nspace
    edge_cross = Ref(0)
    if dir & 1 == 0
        _pidx = ntuple(Val(D)) do i
            v = @inbounds pidx[i]
            if i != _dir
                return v
            end
            v += 1
            @inbounds if v > nspace[i]
                edge_cross[] = 1
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
                edge_cross[] = -1
                return nspace[i]
            end
            return v
        end
    end
    # @assert pidx != _pidx
    return _pidx, edge_cross[]
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
        if step.tidx == param.ntime
            state.nlines += 1
            tidx_next = 1
        else
            tidx_next = step.tidx + 1
        end

        # Next coordinate
        r = floor(Int, rand(T) / param.P_hop)
        if r >= 2 * D
            dirf = 0x1
            dirb_next_new = 0x1
            state.nstraight += 1
        else
            # Forward direction
            dirf = (r + 2) % UInt8
            # Backward direction for the next site.
            # This is the direction along the same axis but reversed.
            dirb_next_new = dirf ⊻ 0x1
            state.nhop += 1
        end
        set_dir(state, (step.tidx, step.pidx...), (dirf, step.dirb))
        pidx_next, edge_cross = get_next(state, step.pidx, dirf)
        state.winding += edge_cross
        dirf_next, dirb_next_orig = get_dir(state, (tidx_next, pidx_next...))
        if dirf_next == 0
            # Empty site, we might want to turn around.
            turn_around = param.forward_turn && rand(T) < param.P_turn
            if !turn_around
                state.χ₀ += 1
            end
            # In either case, the `dirb` for the step
            # needs to be the backward direction we've just created.
            return false, Step{D}(!turn_around, tidx_next, pidx_next, dirb_next_new)
        end
        # We are either turning around (since we hit a line in the middle)
        # or stopping (since we've hit a start).
        # In both cases, we need to save the state for the next site since
        # the next step (if it exists) won't save this.
        set_dir(state, (tidx_next, pidx_next...), (dirf_next, dirb_next_new))
        # Finish if the line we hit didn't have a backward edge (it's a starting point)
        return dirb_next_orig == 0, Step{D}(false, tidx_next, pidx_next, dirb_next_orig)
    else
        # Next time
        tidx_next = step.tidx - 1
        if tidx_next <= 0
            state.nlines -= 1
            tidx_next = param.ntime
        end

        # Next coordinate
        if step.dirb == 1
            state.nstraight -= 1
        else
            state.nhop -= 1
        end
        pidx_next, edge_cross = get_next(state, step.pidx, step.dirb)
        state.winding -= edge_cross
        dirf_next, dirb_next = get_dir(state, (tidx_next, pidx_next...))
        # @assert dirf_next != 0
        set_dir(state, (tidx_next, pidx_next...), (0x0, 0x0))

        turn_around = !param.forward_turn && rand(T) < param.P_turn
        if turn_around
            state.χ₀ += 1
        end
        state.χ₀ += 1
        return (!turn_around && dirb_next == 0,
                Step{D}(turn_around, tidx_next, pidx_next, dirb_next))
    end
end

function sweep!(state::State{T,D}) where {T,D}
    param = state.param
    tidx, pidx... = rand(CartesianIndices((param.ntime, param.nspace...))).I
    dirf, dirb = get_dir(state, (tidx, pidx...))

    step = Step{D}(dirb == 0, tidx, pidx, dirb)

    if step.forward
        state.χ₀ = 1
    else
        state.χ₀ = 0
        # sweep_iterate! backward expect the state to be updated for the site.
        set_dir(state, (tidx, pidx...), (dirf, 0x0))
    end

    while true
        finish, step = sweep_iterate!(state, step)
        if finish
            break
        end
    end

    return (state.nlines,
            state.nstraight * state.param.Estraight + state.nhop * state.param.Ehop,
            state.winding, state.χ₀)
end

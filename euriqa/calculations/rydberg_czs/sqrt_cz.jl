#

using StaticArrays
using LinearAlgebra

@inline function rot(sθ, cθ, sϕ, cϕ)
    sϕ, cϕ = (sϕ, cϕ) .* sθ
    return @SMatrix [cθ complex(sϕ, -cϕ)
                     complex(-sϕ, -cϕ) cθ]
end
@inline function rot_grad(sθ, sϕ, cϕ)
    sϕ, cϕ = (sϕ, cϕ) .* sθ
    return @SMatrix [0 complex(cϕ, sϕ)
                     complex(-cϕ, sϕ) 0]
end

@inline function rot_01(θ, ϕ)
    s, c = sincos(θ / 2)
    sp, cp = sincos(ϕ) .* s
    return @SMatrix [c complex(sp, -cp)
                     complex(-sp, -cp) c]
end

@inline function rot_01_grad(θ, ϕ)
    s = sin(θ / 2)
    sp, cp = sincos(ϕ) .* s
    return @SMatrix [0 complex(cp, sp)
                     complex(-cp, sp) 0]
end

@inline rot_11(θ, ϕ) = rot_01(sqrt(2) * θ, ϕ)
@inline rot_11_grad(θ, ϕ) = rot_01_grad(sqrt(2) * θ, ϕ)


@inline function rydberg_time_wgrad(scθ, scϕ, dt, grads)
    N = length(scθ.s)
    U = @inline ntuple(@inline(i->rot(scθ.s[i], scθ.c[i], scϕ.s[i], scϕ.c[i])), Val(N))
    R0 = @SVector ComplexF64[1, 0]
    Rs = MVector{N+1,typeof(R0)}(undef)
    @inbounds Rs[1] = R0
    R = R0
    tau = 0.0
    @inbounds for i in 1:N
        R = U[i] * R
        Rs[i + 1] = R
        tau += abs2(R[2])
    end
    tau *= dt
    if !isempty(grads)
        e0 = @SVector [1, 0]
        w = @inbounds Rs[N + 1][1] * e0
        @inbounds for j in N:-1:1
            rg = rot_grad(scθ.s[j], scϕ.s[j], scϕ.c[j])
            grads[j] = -2 * dt * real(w' * rg * Rs[j])
            w = muladd.(Rs[j][1], e0, U[j]' * w)
        end
    end
    return tau
end

@inline function infid_sqrtCZ_robust_wgrad(scθ01, scθ11, scϕ, dt, lam, infid_grads)
    N = length(scθ01.s)
    rl0 = @SVector ComplexF64[1, 0]
    right01 = MVector{N,typeof(rl0)}(undef)
    right11 = MVector{N,typeof(rl0)}(undef)
    left01 = MVector{N,typeof(rl0)}(undef)
    left11 = MVector{N,typeof(rl0)}(undef)
    @inbounds begin
        right01[1] = rl0
        right11[1] = rl0
        left01[1] = rl0
        left11[1] = rl0
    end
    @inbounds @simd for i in 1:N - 1
        r01 = rot(scθ01.s[i], scθ01.c[i], scϕ.s[i], scϕ.c[i])
        r11 = rot(scθ11.s[i], scθ11.c[i], scϕ.s[i], scϕ.c[i])
        right01[i + 1] = r01 * right01[i]
        right11[i + 1] = r11 * right11[i]
    end
    @inbounds @simd for i in 1:N - 1
        j = N + 1 - i
        r01 = rot(scθ01.s[j], scθ01.c[j], scϕ.s[j], scϕ.c[j])
        r11 = rot(scθ11.s[j], scθ11.c[j], scϕ.s[j], scϕ.c[j])
        left01[i + 1] = transpose(r01) * left01[i]
        left11[i + 1] = transpose(r11) * left11[i]
    end
    @inbounds begin
        pend = complex(scϕ.c[end], scϕ.s[end])
        pend2 = pend^2

        fid_01 = transpose(left01[2]) * right01[end] * pend
        fid_11 = transpose(left11[2]) * right11[end] * pend2
    end

    A = 1 + 2 * fid_01 - im * fid_11

    infid = 1 - abs2(A) * 0.0625

    if isempty(infid_grads)
        tau01 = rydberg_time_wgrad(scθ01, scϕ, dt, ())
        tau11 = rydberg_time_wgrad(scθ11, scϕ, dt, ())
        d = muladd(-2, tau01, tau11)
        return muladd(lam, abs2(d), infid)
    end

    g01 = MVector{N + 1,ComplexF64}(undef)
    g11 = MVector{N + 1,ComplexF64}(undef)

    @inbounds @simd for i in 1:N
        rg01 = rot_grad(scθ01.s[i], scϕ.s[i], scϕ.c[i])
        rg11 = rot_grad(scθ11.s[i], scϕ.s[i], scϕ.c[i])
        g01[i] = (transpose(left01[N + 1 - i]) * rg01 * right01[i]) * pend
        g11[i] = (transpose(left11[N + 1 - i]) * rg11 * right11[i]) * pend2
    end
    @inbounds begin
        g01[end] = im * fid_01
        g11[end] = 2im * fid_11
    end
    @inbounds @simd for i in 1:N + 1
        infid_grads[i] = -0.125 * real(A' * (2 * g01[i] - im * g11[i]))
    end
    dtau01 = MVector{N, Float64}(undef)
    dtau11 = MVector{N, Float64}(undef)
    tau01 = rydberg_time_wgrad(scθ01, scϕ, dt, dtau01)
    tau11 = rydberg_time_wgrad(scθ11, scϕ, dt, dtau11)

    d = muladd(-2, tau01, tau11)
    # robustness has no Z-phase gradient
    @inbounds @simd for i in 1:N
        infid_grads[i] = muladd(2 * lam * d, muladd(-2, dtau01[i], dtau11[i]),
                                infid_grads[i])
    end
    return muladd(lam, abs2(d), infid)
end

function leak_amp_wgrad(θ_list, pm, @specialize(rot), @specialize(rot_grad),
                        dt, grads)
    N = length(θ_list)
    S2 = @SMatrix [0 0.5;
                   0.5 0]
    U = @inline ntuple(@inline(i->rot(θ_list[i], pm[i])), Val(N))
    R0 = @SVector ComplexF64[1, 0]
    Rs = MVector{N+1,typeof(R0)}(undef)
    @inbounds Rs[1] = R0
    R = R0
    I = 0.0im
    @inbounds for i in 1:N
        R = U[i] * R
        Rs[i + 1] = R
        I = muladd(R[1], R[2], I)
    end
    I *= dt
    if !isempty(grads)
        w = @inbounds S2 * Rs[end]
        @inbounds for j in N:-1:1
            grads[j] = 2 * dt * (transpose(w) * rot_grad(θ_list[j], pm[j]) * Rs[j])
            w = S2 * Rs[j] .+ transpose(U[j]) * w
        end
    end
    return I
end

function darkleak_amp_wgrad(θ_list, pm, @specialize(rot),
                            @specialize(rot_grad), dt, grads)
    N = length(θ_list)
    e1 = @SVector [0, 1]
    U = @inline ntuple(@inline(i->rot(θ_list[i], pm[i])), Val(N))
    R0 = @SVector ComplexF64[1, 0]
    Rs = MVector{N+1,typeof(R0)}(undef)
    @inbounds Rs[1] = R0
    R = R0
    I = 0.0im
    @inbounds for i in 1:N
        R = U[i] * R
        Rs[i + 1] = R
        I += R[2]
    end
    I *= dt
    if !isempty(grads)
        w = @SVector ComplexF64[0, 1]
        @inbounds for j in N:-1:1
            grads[j] = dt * (transpose(w) * rot_grad(θ_list[j], pm[j]) * Rs[j])
            w = e1 .+ transpose(U[j]) * w
        end
    end
    return I
end

function infid_full_wgrad(Ωs, t_gate, ϕs, lam_rob, lam_leak, lam_dark, grads)
    N = length(Ωs)
    @assert length(ϕs) == N + 1
    dt = t_gate / N

    scθ01 = (s = MVector{N,Float64}(undef),
              c = MVector{N,Float64}(undef))
    scθ11 = (s = MVector{N,Float64}(undef),
              c = MVector{N,Float64}(undef))
    @inbounds for i in 1:N
        Ω = Ωs[i]
        θ = Ω * dt
        s01, c01 = sincos(θ / 2)
        s11, c11 = sincos(θ / sqrt(2))
        scθ01.s[i] = s01
        scθ01.c[i] = c01
        scθ11.s[i] = s11
        scθ11.c[i] = c11
    end

    scϕ = (s = MVector{N + 1,Float64}(undef),
           c = MVector{N + 1,Float64}(undef))
    @inbounds for i in 1:N + 1
        scϕ.s[i], scϕ.c[i] = sincos(ϕs[i])
    end

    J = infid_sqrtCZ_robust_wgrad(scθ01, scθ11, scϕ, dt, lam_rob, grads)
    θs = Ωs .* dt
    pm = ϕs
    if isempty(grads)
        dI01 = ()
        dI11 = ()
        dId11 = ()
    else
        dI01 = MVector{N, ComplexF64}(undef)
        dI11 = MVector{N, ComplexF64}(undef)
        dId11 = MVector{N, ComplexF64}(undef)
    end
    I01 = leak_amp_wgrad(θs, pm, rot_01, rot_01_grad, dt, dI01)
    I11 = leak_amp_wgrad(θs, pm, rot_11, rot_11_grad, dt, dI11)
    Id11 = darkleak_amp_wgrad(θs, pm, rot_11, rot_11_grad, dt, dId11)
    J = muladd(lam_leak, muladd(2, abs2(I01), abs2(I11)), muladd(lam_dark, abs2(Id11), J))
    if !isempty(grads)
        @inbounds for i in 1:N
            grads[i] += lam_leak * (4 * real(I01' * dI01[i]) + 2 * real(I11' * dI11[i])) + lam_dark * 2 * real(Id11' * dId11[i])
        end
    end
    return J
end

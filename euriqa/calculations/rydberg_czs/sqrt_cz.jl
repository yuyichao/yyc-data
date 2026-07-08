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

const SM2{T} = SMatrix{2, 2, T, 4}
struct InfidFullData{N,has_grad,N2,TGR}
    U01::MVector{N,SM2{ComplexF64}}
    U11::MVector{N,SM2{ComplexF64}}
    R01::MVector{N2,SVector{2,ComplexF64}}
    R11::MVector{N2,SVector{2,ComplexF64}}
    gR01::TGR
    gR11::TGR
    pend::ComplexF64
    pend2::ComplexF64
    dt::Float64
    @inline function InfidFullData(Ωs, ϕs, t_gate, ::Val{has_grad}) where has_grad
        N = length(Ωs)
        N2 = N + 1
        @assert length(ϕs) == N2
        dt = t_gate / N
        U01 = MVector{N,SM2{ComplexF64}}(undef)
        U11 = MVector{N,SM2{ComplexF64}}(undef)
        R01 = MVector{N2,SVector{2,ComplexF64}}(undef)
        R11 = MVector{N2,SVector{2,ComplexF64}}(undef)
        if has_grad
            gR01 = MVector{N,SVector{2,ComplexF64}}(undef)
            gR11 = MVector{N,SVector{2,ComplexF64}}(undef)
        else
            gR01 = ()
            gR11 = ()
        end
        TGR = typeof(gR01)

        r01 = @SVector ComplexF64[1, 0]
        r11 = @SVector ComplexF64[1, 0]
        @inbounds begin
            R01[1] = r01
            R11[1] = r11
        end
        @inbounds for i in 1:N
            Ω = Ωs[i]
            θ = Ω * dt
            s01, c01 = sincos(θ / 2)
            s11, c11 = sincos(θ / sqrt(2))
            sϕ, cϕ = sincos(ϕs[i])
            u01 = rot(s01, c01, sϕ, cϕ)
            u11 = rot(s11, c11, sϕ, cϕ)
            if has_grad
                g01 = rot_grad(s01, sϕ, cϕ)
                g11 = rot_grad(s11, sϕ, cϕ)
                gR01[i] = g01 * R01[i]
                gR11[i] = g11 * R11[i]
            end
            U01[i] = u01
            U11[i] = u11

            r01 = u01 * r01
            r11 = u11 * r11
            R01[i + 1] = r01
            R11[i + 1] = r11
        end
        pend = @inbounds cis(ϕs[N + 1])
        pend2 = pend^2
        return new{N,has_grad,N2,TGR}(U01, U11, R01, R11, gR01, gR11, pend, pend2, dt)
    end
end

@inline function rydberg_time_wgrad(Us, Rs, gRs, dt, grads)
    N = length(Us)
    tau = 0.0
    @inbounds for i in 2:N + 1
        tau += abs2(Rs[i][2])
    end
    tau *= dt
    if grads !== ()
        e0 = @SVector [1, 0]
        w = @inbounds Rs[N + 1][1] * e0
        @inbounds for j in N:-1:1
            grads[j] = -2 * dt * real(w' * gRs[j])
            w = muladd.(Rs[j][1], e0, Us[j]' * w)
        end
    end
    return tau
end

@inline function infid_sqrtCZ_robust_wgrad(d::InfidFullData{N,has_grad}, lam, infid_grads) where {N,has_grad}
    dt = d.dt
    l0 = @SVector ComplexF64[1, 0]
    left01 = MVector{N,typeof(l0)}(undef)
    left11 = MVector{N,typeof(l0)}(undef)
    @inbounds begin
        left01[1] = l0
        left11[1] = l0
    end
    @inbounds @simd for i in 1:N - 1
        j = N + 1 - i
        left01[i + 1] = transpose(d.U01[j]) * left01[i]
        left11[i + 1] = transpose(d.U11[j]) * left11[i]
    end
    @inbounds begin
        fid_01 = transpose(left01[2]) * d.R01[N] * d.pend
        fid_11 = transpose(left11[2]) * d.R11[N] * d.pend2
    end

    A = 1 + 2 * fid_01 - im * fid_11

    infid = 1 - abs2(A) * 0.0625

    if !has_grad
        tau01 = rydberg_time_wgrad(d.U01, d.R01, d.gR01, dt, ())
        tau11 = rydberg_time_wgrad(d.U11, d.R11, d.gR11, dt, ())
        d = muladd(-2, tau01, tau11)
        return muladd(lam, abs2(d), infid)
    end

    g01 = MVector{N + 1,ComplexF64}(undef)
    g11 = MVector{N + 1,ComplexF64}(undef)

    @inbounds @simd for i in 1:N
        g01[i] = (transpose(left01[N + 1 - i]) * d.gR01[i]) * d.pend
        g11[i] = (transpose(left11[N + 1 - i]) * d.gR11[i]) * d.pend2
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
    tau01 = rydberg_time_wgrad(d.U01, d.R01, d.gR01, dt, dtau01)
    tau11 = rydberg_time_wgrad(d.U11, d.R11, d.gR11, dt, dtau11)

    d = muladd(-2, tau01, tau11)
    # robustness has no Z-phase gradient
    @inbounds @simd for i in 1:N
        infid_grads[i] = muladd(2 * lam * d, muladd(-2, dtau01[i], dtau11[i]),
                                infid_grads[i])
    end
    return muladd(lam, abs2(d), infid)
end

@inline function leak_amp_wgrad(Us, Rs, gRs, dt, grads)
    N = length(Us)
    I = 0.0im
    @inbounds for i in 2:N + 1
        R = Rs[i]
        I = muladd(R[1], R[2], I)
    end
    I *= dt
    if grads !== ()
        S2 = @SMatrix [0 0.5
                       0.5 0]
        w = @inbounds S2 * Rs[end]
        @inbounds for j in N:-1:1
            grads[j] = 2 * dt * (transpose(w) * gRs[j])
            w = S2 * Rs[j] .+ transpose(Us[j]) * w
        end
    end
    return I
end

@inline function darkleak_amp_wgrad(Us, Rs, gRs, dt, grads)
    N = length(Us)
    I = 0.0im
    @inbounds for i in 2:N + 1
        I += Rs[i][2]
    end
    I *= dt
    if grads !== ()
        e1 = @SVector [0, 1]
        w = @SVector ComplexF64[0, 1]
        @inbounds for j in N:-1:1
            grads[j] = dt * (transpose(w) * gRs[j])
            w = e1 .+ transpose(Us[j]) * w
        end
    end
    return I
end

function infid_full_wgrad(Ωs, t_gate, ϕs, lam_rob, lam_leak, lam_dark, grads)
    N = length(Ωs)
    @assert length(ϕs) == N + 1
    d = InfidFullData(Ωs, ϕs, t_gate, Val(!isempty(grads)))
    J = infid_sqrtCZ_robust_wgrad(d, lam_rob, grads)

    dt = t_gate / N

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
    I01 = leak_amp_wgrad(d.U01, d.R01, d.gR01, dt, dI01)
    I11 = leak_amp_wgrad(d.U11, d.R11, d.gR11, dt, dI11)
    Id11 = darkleak_amp_wgrad(d.U11, d.R11, d.gR11, dt, dId11)
    J = muladd(lam_leak, muladd(2, abs2(I01), abs2(I11)), muladd(lam_dark, abs2(Id11), J))
    if !isempty(grads)
        @inbounds for i in 1:N
            grads[i] += lam_leak * (4 * real(I01' * dI01[i]) + 2 * real(I11' * dI11[i])) + lam_dark * 2 * real(Id11' * dId11[i])
        end
    end
    return J
end

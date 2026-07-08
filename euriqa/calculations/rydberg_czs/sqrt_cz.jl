#

using StaticArrays
using Static
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
            s01, c01 = @fastmath sincos(θ / 2)
            s11, c11 = @fastmath sincos(θ / sqrt(2))
            sϕ, cϕ = @fastmath sincos(ϕs[i])
            u01 = rot(s01, c01, sϕ, cϕ)
            u11 = rot(s11, c11, sϕ, cϕ)
            if has_grad
                gR01[i] = rot_grad(s01, sϕ, cϕ) * r01
                gR11[i] = rot_grad(s11, sϕ, cϕ) * r11
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

@inline function rydberg_time_wgrad(Us, Rs, gRs, dt, grads, gα, gβ)
    N = length(Us)
    tau = 0.0
    @inbounds for i in 2:N + 1
        tau += abs2(Rs[i][2])
    end
    tau *= dt
    if grads !== ()
        e0 = @SVector [1, 0]
        w = @inbounds Rs[N + 1][1] * e0
        if gα == 0
            @inbounds for j in N:-1:1
                grads[j] = -2 * dt * gβ * real(w' * gRs[j])
                w = muladd.(Rs[j][1], e0, Us[j]' * w)
            end
        else
            @inbounds for j in N:-1:1
                grads[j] = muladd(-2 * dt * gβ, real(w' * gRs[j]), grads[j] * gα)
                w = muladd.(Rs[j][1], e0, Us[j]' * w)
            end
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
    @inbounds for i in 1:N - 1
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
        if lam > 0
            tau01 = rydberg_time_wgrad(d.U01, d.R01, d.gR01, dt,
                                       (), static(false), static(false))
            tau11 = rydberg_time_wgrad(d.U11, d.R11, d.gR11, dt,
                                       (), static(false), static(false))
            d = muladd(-2, tau01, tau11)
            infid = muladd(lam, abs2(d), infid)
        end
        return infid
    end

    if lam > 0
        tau01 = rydberg_time_wgrad(d.U01, d.R01, d.gR01, dt,
                                   infid_grads, static(false), static(-2))
        tau11 = rydberg_time_wgrad(d.U11, d.R11, d.gR11, dt,
                                   infid_grads, static(true), static(true))
        tau = muladd(-2, tau01, tau11)

        @inbounds for i in 1:N
            g01 = (transpose(left01[N + 1 - i]) * d.gR01[i]) * d.pend
            g11 = (transpose(left11[N + 1 - i]) * d.gR11[i]) * d.pend2
            g11 = complex(imag(g11), -real(g11))
            dtau = infid_grads[i]
            infid_grads[i] = muladd(2 * lam * tau, dtau,
                                    -0.125 * real(A' * muladd(2, g01, g11)))
        end
        infid = muladd(lam, abs2(tau), infid)
    else
        @inbounds for i in 1:N
            g01 = (transpose(left01[N + 1 - i]) * d.gR01[i]) * d.pend
            g11 = (transpose(left11[N + 1 - i]) * d.gR11[i]) * d.pend2
            g11 = complex(imag(g11), -real(g11))
            infid_grads[i] = -0.125 * real(A' * muladd(2, g01, g11))
        end
    end
    @inbounds begin
        g01 = im * fid_01
        g11 = 2im * fid_11
        infid_grads[N + 1] = -0.125 * real(A' * (2 * g01 - im * g11))
    end
    return infid
end

@inline function leak_amp_wgrad(Us, Rs, gRs, dt, grads, gβ, ::Val{has_grad}) where has_grad
    N = length(Us)
    I = 0.0im
    @inbounds for i in 2:N + 1
        R = Rs[i]
        I = muladd(R[1], R[2], I)
    end
    I *= dt
    if has_grad
        function S2(R)
            @SVector [R[2], R[1]]
        end
        w = @inbounds S2(Rs[end]) .* 0.5
        @inbounds for j in N:-1:1
            g = transpose(w) * gRs[j]
            grads[j] = muladd(4 * gβ * dt, real(I' * g), grads[j])
            w = muladd.(S2(Rs[j]), 0.5, transpose(Us[j]) * w)
        end
    end
    return abs2(I)
end

@inline function darkleak_amp_wgrad(Us, Rs, gRs, dt, grads, gβ, ::Val{has_grad}) where has_grad
    N = length(Us)
    I = 0.0im
    @inbounds for i in 2:N + 1
        I += Rs[i][2]
    end
    I *= dt
    if grads !== ()
        e1 = @SVector [false, true]
        w = @SVector ComplexF64[0, 1]
        @inbounds for j in N:-1:1
            g = transpose(w) * gRs[j]
            grads[j] = muladd(2 * gβ * dt, real(I' * g), grads[j])
            w = e1 .+ transpose(Us[j]) * w
        end
    end
    return abs2(I)
end

@inline function infid_full_wgrad(d::InfidFullData{N,has_grad},
                                  lam_rob, lam_leak, lam_dark, grads) where {N,has_grad}
    dt = d.dt
    J = infid_sqrtCZ_robust_wgrad(d, lam_rob, grads)
    if lam_leak > 0
        I01 = leak_amp_wgrad(d.U01, d.R01, d.gR01, dt, grads, 2 * lam_leak, Val(has_grad))
        I11 = leak_amp_wgrad(d.U11, d.R11, d.gR11, dt, grads, lam_leak, Val(has_grad))
        J = muladd(lam_leak, muladd(2, I01, I11), J)
    end
    if lam_dark > 0
        Id11 = darkleak_amp_wgrad(d.U11, d.R11, d.gR11, dt, grads, lam_dark, Val(has_grad))
        J = muladd(lam_dark, Id11, J)
    end
    return J
end

function infid_full_wgrad(Ωs, t_gate, ϕs, lam_rob, lam_leak, lam_dark, grads)
    N = length(Ωs)
    @assert length(ϕs) == N + 1
    if isempty(grads)
        return infid_full_wgrad(InfidFullData(Ωs, ϕs, t_gate, Val(false)),
                                lam_rob, lam_leak, lam_dark, ())
    else
        return infid_full_wgrad(InfidFullData(Ωs, ϕs, t_gate, Val(true)),
                                lam_rob, lam_leak, lam_dark, grads)
    end
end

#

using StaticArrays
using Static
using LinearAlgebra
using NLopt
using MSSim: Optimizers as Opts

@inline function rot(sθ, cθ, sϕ, cϕ)
    sϕ, cϕ = (sϕ, cϕ) .* sθ
    return @SMatrix [cθ complex(sϕ, -cϕ)
                     complex(-sϕ, -cϕ) cθ]
end
@inline function rot_grad(sθ, sϕ, cϕ, r)
    sϕ, cϕ = (sϕ, cϕ) .* sθ
    return @inbounds @SVector [complex(cϕ, sϕ) * r[2],
                               complex(-cϕ, sϕ) * r[1]]
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
    @inline function InfidFullData(Ωs, ϕs, dt, ::Val{has_grad}) where has_grad
        N = length(Ωs)
        N2 = N + 1
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
                gR01[i] = rot_grad(s01, sϕ, cϕ, r01)
                gR11[i] = rot_grad(s11, sϕ, cϕ, r11)
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
            @inbounds @SVector [R[2], R[1]]
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

@inline function infid_full_wgrad(Ωs, dt, ϕs, lam_rob, lam_leak, lam_dark, grads)
    N = length(Ωs)
    @assert length(ϕs) == N + 1
    if isempty(grads)
        return infid_full_wgrad(InfidFullData(Ωs, ϕs, dt, Val(false)),
                                lam_rob, lam_leak, lam_dark, ())
    else
        return infid_full_wgrad(InfidFullData(Ωs, ϕs, dt, Val(true)),
                                lam_rob, lam_leak, lam_dark, grads)
    end
end

function spline_matrix(nseg, nsample)
    MD = zeros(nseg + 1, nseg + 1)
    @inbounds for i in 1:nseg
        MD[i, i] = 4
        MD[i, i + 1] = 1
        MD[i + 1, i] = 1
    end
    @inbounds MD[1, 1] = 2
    @inbounds MD[nseg + 1, nseg + 1] = 2
    M3 = zeros(nseg + 1, nseg + 1)
    @inbounds for i in 2:nseg
        M3[i, i + 1] = 3
        M3[i, i - 1] = -3
    end
    @inbounds M3[1, 2] = M3[nseg + 1, nseg + 1] = 3
    @inbounds M3[1, 1] = M3[nseg + 1, nseg] = -3
    SD = MD \ M3
    Mabcd = zeros(4 * nseg, nseg + 1)
    @inbounds for i in 1:nseg
        offset = (i - 1) * 4
        # aᵢ = yᵢ
        Mabcd[offset + 1, i] = 1
        # bᵢ = Dᵢ
        Mabcd[offset + 2, :] .= @view(SD[i, :])
        # cᵢ = 3(yᵢ₊₁ - yᵢ) - 2Dᵢ - Dᵢ₊₁
        Mabcd[offset + 3, :] .= muladd.(-2, @view(SD[i, :]), .-@view(SD[i + 1, :]))
        Mabcd[offset + 3, i] -= 3
        Mabcd[offset + 3, i + 1] += 3
        # dᵢ = 2(yᵢ - yᵢ₊₁) + Dᵢ + Dᵢ₊₁
        Mabcd[offset + 4, :] .= @view(SD[i, :]) .+ @view(SD[i + 1, :])
        Mabcd[offset + 4, i] += 2
        Mabcd[offset + 4, i + 1] -= 2
    end
    M = zeros(nseg * nsample + 1, nseg + 1)
    @inbounds for i in 1:nseg
        offset_abcd = (i - 1) * 4
        offset = (i - 1) * nsample
        for j in 1:(i == nseg ? nsample + 1 : nsample)
            x = (j - 1) / nsample
            a = @view(Mabcd[offset_abcd + 1, :])
            b = @view(Mabcd[offset_abcd + 2, :])
            c = @view(Mabcd[offset_abcd + 3, :])
            d = @view(Mabcd[offset_abcd + 4, :])
            M[offset + j, :] .= a .+ x .* (b .+ x .* (c .+ x .* d))
        end
    end
    return M
end

function fm_spline_matrix(ttotal, nseg, nsample)
    dt = ttotal / (nseg * nsample)
    Ms = spline_matrix(nseg, nsample)
    M = zeros(nseg * nsample + 1, nseg + 1)
    @inbounds for i in 1:nseg + 1
        v = 0.0
        for j in 2:nseg * nsample + 1
            v += Ms[j - 1, i] * dt
            M[j, i] = v
        end
    end
    return M
end

@inline function infid_full_fm(Ωs, dt, ωsϕ, lam_rob, lam_leak, lam_dark, grads)
    N = length(Ωs)
    @assert length(ωsϕ) == N
    ϕs = MVector{N + 1,Float64}(undef)
    ϕgrads = MVector{N + 1,Float64}(undef)
    @inbounds ϕs[1] = 0
    ϕ = 0.0
    @inbounds for i in 1:N - 1
        ϕ = muladd(dt, ωsϕ[i], ϕ)
        ϕs[i + 1] = ϕ
    end
    @inbounds ϕs[N + 1] = ωsϕ[N]
    infid = infid_full_wgrad(Ωs, dt, ϕs, lam_rob, lam_leak, lam_dark, ϕgrads)
    @inbounds if !isempty(grads)
        g = 0.0
        for i in N - 1:-1:1
            g = muladd(ϕgrads[i + 1], dt, g)
            grads[i] = g
        end
        grads[N] = ϕgrads[N + 1]
    end
    return infid
end

function get_infid_full_cb(Ωs, t_gate, lam_rob, lam_leak, lam_dark)
    dt = t_gate / length(Ωs)
    return function infid_full_cb(ϕs, grads)
        return infid_full_wgrad(Ωs, dt, ϕs, lam_rob, lam_leak, lam_dark, grads)
    end
end

function get_infid_full_fm_cb(Ωs, t_gate, lam_rob, lam_leak, lam_dark)
    dt = t_gate / length(Ωs)
    return function infid_full_cb(ωsϕ, grads)
        return infid_full_fm(Ωs, dt, ωsϕ, lam_rob, lam_leak, lam_dark, grads)
    end
end

fm_to_phase(ωsϕ, t_gate) = [0; cumsum(@view(ωsϕ[1:end - 1])) .* (t_gate / length(ωsϕ));
                             ωsϕ[end]]

struct Opt{CB}
    opt::NLopt.Opt
    pre_opt::NLopt.Opt
    tracker::Opts.NLVarTracker
    args_buff::Vector{Float64}
    cb::CB
    t_gate::Float64
end

function Opt(Ωs, t_gate, lam_rob, lam_leak, lam_dark;
             algorithm=:LD_CCSAQ, maxeval_pre=1000,
             maxtime=3, xtol=1e-7, minω=-2π * 10, maxω=2π * 10)
    N = length(Ωs)
    tracker = Opts.NLVarTracker(N)
    opt = NLopt.Opt(algorithm, N)
    pre_opt = NLopt.Opt(algorithm, N)

    for i in 1:N - 1
        Opts.set_bound!(tracker, i, minω, maxω)
    end
    Opts.set_bound!(tracker, N, -4π, 4π)

    NLopt.maxeval!(pre_opt, maxeval_pre)
    NLopt.xtol_rel!(pre_opt, xtol * 10)
    NLopt.lower_bounds!(pre_opt, Opts.lower_bounds(tracker))
    NLopt.upper_bounds!(pre_opt, Opts.upper_bounds(tracker))

    NLopt.maxtime!(opt, maxtime)
    NLopt.xtol_rel!(opt, xtol)
    NLopt.lower_bounds!(opt, Opts.lower_bounds(tracker))
    NLopt.upper_bounds!(opt, Opts.upper_bounds(tracker))
    cb = get_infid_full_fm_cb(Ωs, t_gate, static(lam_rob),
                              static(lam_leak), static(lam_dark))
    NLopt.min_objective!(pre_opt, cb)
    NLopt.min_objective!(opt, cb)
    return Opt{typeof(cb)}(opt, pre_opt, tracker, Vector{Float64}(undef, N), cb, t_gate)
end

function Opt(; Ω, num_slices, t_gate, lam_rob, lam_leak, lam_dark, kws...)
    Ωs = @SVector(fill(Ω, num_slices))
    return Opt(Ωs, t_gate, lam_rob, lam_leak, lam_dark; kws...)
end

fm_to_phase(opt::Opt, ωsϕ) = fm_to_phase(ωsϕ, opt.t_gate)

function opt_one!(opt::Opt, pre_threshold)
    objval, args, ret = NLopt.optimize!(opt.pre_opt,
                                        Opts.init_vars!(opt.tracker, opt.args_buff))
    if getfield(NLopt, ret)::NLopt.Result < 0
        return
    end
    if objval > pre_threshold
        return objval, args
    end
    objval, args, ret = NLopt.optimize!(opt.opt, args)
    if getfield(NLopt, ret)::NLopt.Result < 0
        return
    end
    return objval, args
end

function opt_n!(opt::Opt, n; verbose=true, pre_threshold=0.005,
                best_obj=1.0, best_args=similar(opt.args_buff))
    for i in 1:n
        if verbose
            opt_res = @time opt_one!(opt, pre_threshold)
        else
            opt_res = opt_one!(opt, pre_threshold)
            if i % 20 == 0
                println("Round $i done")
            end
        end
        if opt_res === nothing
            continue
        end
        obj, args = opt_res
        if obj < best_obj
            @show obj
            best_obj = obj
            best_args .= args
        end
    end
    return best_obj, best_args
end

# const opt = get_opt(Ω=2π * 3, num_slices=100, t_gate=0.6,
#                     lam_rob=0, lam_leak=1, lam_dark=1)

# Ωs = @SVector [0.8120548277768912, 0.6600756672118772, 0.8834624513620218, 0.3784616113155762, 0.718592585456476, 0.2157486028591974, 0.4759595586914189, 0.8853953650900385, 0.38569591634297806, 0.5425073635198567, 0.4547036110530509, 0.36524492058833724, 0.9336343529283576, 0.5417206805604657, 0.4264683528030687, 0.7436124287466458, 0.9003200753748883, 0.6194226364448447, 0.8905048715742345, 0.41242634095552677, 0.8327270563601665, 0.7028116134210601, 0.10351802031716417, 0.6708798796477251, 0.4878743133421236, 0.1326162534712254, 0.2905364957596993, 0.7116240123204844, 0.9291825808455354, 0.032078975375203655, 0.2286254558149481, 0.25430294718459723]
# ϕs = @SVector [0.9778373273164378, 0.06357512151726186, 0.9558640651034589, 0.6471022152045088, 0.12347033934226914, 0.022044412236653654, 0.9058383482270784, 0.7628439979061399, 0.2447784047455075, 0.4631491990304336, 0.17609449014012069, 0.05016441975120223, 0.8191325108006975, 0.01312912412594236, 0.7116411798906169, 0.18799971307390095, 0.47137853828664955, 0.8909624415832684, 0.40548870824966243, 0.849193421455286, 0.6495156770548955, 0.21323705504676804, 0.21681191098550467, 0.8758800971626411, 0.5416163718880839, 0.6805477197934157, 0.9305298996463732, 0.6558294963939725, 0.08138627788491504, 0.8996129066308572, 0.654624115568562, 0.81096160781888, 0.8738053956553339]
# infid_grads = MVector{33, Float64}(undef)
# infid = infid_full_wgrad(Ωs, 2.3 / 32, ϕs, 0.1, 0.2, 0.3, infid_grads)

# using Test
# @test infid ≈ 1.03535111363171 atol=1e-7
# @test infid_grads ≈ [-0.048287795862711896, 0.041540304006545015, -0.04661859027492562, -0.004924112095923546, 0.03598744856862636, 0.012695891004438475, -0.01891394663410687, -0.02053237966215333, 0.011991944547006608, 0.004681818120394754, 0.015825816648202536, 0.015916041484574827, -0.021391491757655527, 0.022831758063062018, -0.005359784647926516, 0.019496742256813304, 0.004678660699337808, -0.014168838597348245, 0.007398012778739823, -0.007668962712466331, -0.005630874839035227, 0.011599665880002673, 0.0015781701677002095, -0.011295049147817322, -0.0006445651991900822, -0.0009705514322887574, -0.00483519052725263, -0.00415453267505124, 0.01338467792380522, -0.0004331684937992036, -0.0012539987131235059, -0.002523118876473579, 0.10274642289173286] atol=1e-7

# using BenchmarkTools
# @btime infid_full_wgrad($(Ωs), $(2.3 / 32), $(ϕs), $(0.1), $(0.2), $(0.3), ())
# @btime infid_full_wgrad($(Ωs), $(2.3 / 32), $(ϕs), $(0.1), $(0.2), $(0.3), $(infid_grads))

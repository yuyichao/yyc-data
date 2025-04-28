#/usr/bin/julia

include("grape.jl")

using JuMP
using NLopt
using StaticArrays

const OPType = SMatrix{4,4,ComplexF64,16}

mutable struct Drive <: AbstractOP{OPType,2}
    θ::Float64
    ϕ::Float64
    Drive() = new()
end
function set_params(dri::Drive, params)
    dri.θ = params[1]
    dri.ϕ = params[2]
    return
end
function compute(dri::Drive, grad)
    st, ct = sincos(dri.θ)
    st2, ct2 = sincos(√2 * dri.θ)

    sp, cp = sincos(dri.ϕ)

    diag1 = ct
    diag1_θ = -st
    diag1_ϕ = 0
    off1 = st * complex(sp, -cp)
    off1_θ = ct * complex(sp, -cp)
    off1_ϕ = st * complex(cp, sp)
    coff1 = st * complex(-sp, -cp)
    coff1_θ = ct * complex(-sp, -cp)
    coff1_ϕ = st * complex(-cp, sp)

    diag2 = ct2
    diag2_θ = -√2 * st2
    diag2_ϕ = 0
    off2 = st2 * complex(sp, -cp)
    off2_θ = √2 * ct2 * complex(sp, -cp)
    off2_ϕ = st2 * complex(cp, sp)
    coff2 = st2 * complex(-sp, -cp)
    coff2_θ = √2 * ct2 * complex(-sp, -cp)
    coff2_ϕ = st2 * complex(-cp, sp)

    grad[1] = @SMatrix[diag1_θ  off1_θ   0        0
                       coff1_θ  diag1_θ  0        0
                       0        0        diag2_θ  off2_θ
                       0        0        coff2_θ  diag2_θ]
    grad[2] = @SMatrix[diag1_ϕ  off1_ϕ   0        0
                       coff1_ϕ  diag1_ϕ  0        0
                       0        0        diag2_ϕ  off2_ϕ
                       0        0        coff2_ϕ  diag2_ϕ]
    return @SMatrix[diag1  off1   0      0
                    coff1  diag1  0      0
                    0      0      diag2  off2
                    0      0      coff2  diag2]
end

mutable struct PhaseZ <: AbstractOP{OPType,1}
    ϕ::Float64
    PhaseZ() = new()
end
function set_params(p::PhaseZ, params)
    p.ϕ = params[1]
    return
end
function compute(p::PhaseZ, grad)
    s, c = sincos(p.ϕ)
    s2 = 2 * s * c
    c2 = 2 * c * c - 1

    diag1 = complex(c, s)
    diag1_ϕ = complex(-s, c)

    diag2 = complex(c2, s2)
    diag2_ϕ = 2 * complex(-s2, c2)

    grad[1] = @SMatrix[diag1_ϕ  0  0        0
                       0        0  0        0
                       0        0  diag2_ϕ  0
                       0        0  0        0]
    return @SMatrix[diag1  0  0      0
                    0      1  0      0
                    0      0  diag2  0
                    0      0  0      1]
end

@inline function convert_res(op)
    return abs2(op[1, 1] - 1) * 2 + abs2(op[3, 3] + 1)
end
@inline function convert_grad(op, op_grad)
    op1 = op[1, 1]
    op3 = op[3, 3]
    op1_grad = op_grad[1, 1]
    op3_grad = op_grad[3, 3]

    return (4 * ((real(op1) - 1) * real(op1_grad) + imag(op1) * imag(op1_grad)) +
        2 * ((real(op3) + 1) * real(op3_grad) + imag(op3) * imag(op3_grad)))
end

@inline function convert_res_grads(op, op_grads, grads)
    grads .= convert_grad.(Ref(op), op_grads)
    return convert_res(op)
end

mutable struct PMPulseSeq{N,S,P,OB,RB}
    const s::S # Sequence
    const params::P # Input parameters
    res::Float64 # Output result
    const grads::P # Output gradients

    const op_buff::OB
    const res_buff::RB

    function PMPulseSeq{N}() where N
        ops = ntuple(Val(N + 1)) do i
            return i == 1 ? PhaseZ() : Drive()
        end
        s = Sequence{OPType}(ops)
        params = MVector{N + 2,Float64}(undef)
        grads = MVector{N + 2,Float64}(undef)
        op_buff = Vector{OPType}(undef,2N + 1)
        res_buff = MVector{2N + 1,Float64}(undef)
        return new{N,typeof(s),typeof(params),typeof(op_buff),typeof(res_buff)}(
            s, params, NaN, grads, op_buff, res_buff)
    end
end

@inline function _set_params_buff!(res_buff, params, N)
    @inbounds begin
        res_buff[1] = params[1]
        angle = params[2] / N
        for i in 1:N
            res_buff[2 * i] = angle
            res_buff[2 * i + 1] = params[i + 2]
        end
    end
end

function update_params!(ps::PMPulseSeq{N}, params) where N
    if !isnan(ps.res) && all(ps.params .== params)
        return
    end
    @assert length(params) == N + 2
    res_buff = ps.res_buff

    _set_params_buff!(res_buff, params, N)
    set_params(ps.s, res_buff)
    op = compute(ps.s, ps.op_buff)
    res = convert_res_grads(op, ps.op_buff, res_buff)

    @inbounds begin
        grads = ps.grads
        grads[1] = res_buff[1]
        angle_grad = 0.0
        for i in 1:N
            angle_grad += res_buff[2 * i]
            grads[i + 2] = res_buff[2 * i + 1]
        end
        grads[2] = angle_grad / N
        ps.res = res
        ps.params .= params
    end
    return
end

function optimize_pulse!(ps::PMPulseSeq{N}, init_params) where N
    m = Model(NLopt.Optimizer)
    set_attribute(m, "algorithm", :LD_SLSQP)
    function res_func(params...)
        @assert length(params) == N + 2
        update_params!(ps, params)
        return ps.res
    end
    function grad_func(g, params...)
        @assert length(params) == N + 2
        update_params!(ps, params)
        g .= ps.grads
        return
    end
    register(m, :fidelity, N + 2, res_func, grad_func, autodiff=false)
    @variable(m, global_z, start=init_params[1])
    @variable(m, total_angle, start=init_params[2])
    @variable(m, phases[i=1:N], start=init_params[i + 2])
    @NLobjective(m, Min, fidelity(global_z, total_angle, phases...))
    JuMP.optimize!(m)
    params = [value(global_z); value(total_angle); value.(phases)]
    res = res_func(params...)
    # @show res
    return res, params
end

# [-4.082623222999641, 4.541397221247234, -0.03245830600737785, -0.8124489994546822, 0.47202474200536804, -0.5500803261909537, -0.6177047316919362, 0.5411791619982704, -0.2557655640738533, -0.38314193680002606, 0.3118226680374636, 0.5029169562968007, -0.19251679080546885, -0.6830286871769535, -0.9731886766315152, 0.44235641455254676, 0.22454690222112847, 0.6042372041827836, 0.5738136841262153, -0.6144887991414903, 0.3956674052845687, -0.4650998175643587, 0.7887892642909274, 0.11058986980661852, 0.3804075261097273, 1.1650585873558115, 0.1872911377975865, 0.9567278306607799, -0.4947561155344319, 0.6894623797962316, 0.09709192320588979, 0.7557011190748159, 0.8259314495280358, 0.8185719680309294, 0.9349236332175126, 1.3108418952354415, 0.8231074271575507, 1.3398388812930984, 1.4380091588537953, 1.3046589025558661, 1.0254758451006798, 0.5774038163683318, 0.8788491987152173, -0.07075782181316924, -0.5112282797272782, 0.6423404924404528, 0.027051359134414, 0.6404011131017019, -0.13714724496857644, -0.6502383866425321, 0.6573836386387235, -0.1668395082070339, 0.3338638845614989, -0.6685430326154859, -0.21722996509706732, 0.06410340375354018, -1.1336402705064414, -1.1331811610610962, -1.1679161618322564, -1.1513218730649466, -0.6240970802227468, -1.1024089561607755, -1.1884290426757993, 0.09776880469789495, -0.3297167756824757, 0.09422865956652068, 0.43128269240808786, 0.27626678345552685, -0.8704769881183879, -1.0788257414112599, 0.06424796464472106, -0.3954332707109115, 0.4754886348364921, -0.2384744955069249, -0.7489516745171744, -0.8813498815632604, -1.164959134556153, 0.577457335488964, -1.2008233920845337, -0.3581232596135072, -0.9058520582238452, -0.7263428645537944, 0.6411906806732239, -1.0437935842604629, -0.7276694585207621, 0.2210376247854085, 0.07694243354081644, 0.44670224404186903, -0.6675923734961656, -0.7336651023078886, -0.30899172164173766, -0.6278027558000905, -0.6015399266994913, 0.19024846578435994, 0.8059652093292294, 0.9677410073384545, 0.26660466658662235, 0.5748930920838033, 0.5499026608622101, -0.23602798130692976, 1.030334485675284, 1.1253630804302641]

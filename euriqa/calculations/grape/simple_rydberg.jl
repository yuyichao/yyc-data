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

mutable struct DetDrive <: AbstractOP{OPType,3}
    θ::Float64
    ϕ::Float64
    α::Float64
    DetDrive() = new()
end
function set_params(dri::DetDrive, params)
    dri.θ = params[1] * 2
    dri.ϕ = params[2]
    dri.α = params[3]
    return
end
@inline function detdrive_terms(θ, ϕ, α, s)
    θ = θ * s

    θ′ = sqrt(abs2(θ) + abs2(α))
    if θ′ == 0
        θ′_θ = 0.0
        θ′_α = 0.0
    else
        θ′_θ = θ / θ′
        θ′_α = α / θ′
    end

    st, ct = sincos(θ′ / 2)
    sct = sinc(θ′ / (2π))
    cct = cosc(θ′ / (2π)) / (2π)

    sa, ca = sincos(α / 2)
    sp, cp = sincos(α / 2 + ϕ)

    dgr = ct * ca + α / 2 * sct * sa
    dgr_θ′ = -st * ca / 2 + α / 2 * cct * sa
    dgr_θ = dgr_θ′ * θ′_θ * s * 2
    dgr_α = dgr_θ′ * θ′_α - ct * sa / 2 + 1 / 2 * sct * sa + α / 4 * sct * ca

    dgi = ct * sa - α / 2 * sct * ca
    dgi_θ′ = -st * sa / 2 - α / 2 * cct * ca
    dgi_θ = dgi_θ′ * θ′_θ * s * 2
    dgi_α = dgi_θ′ * θ′_α + ct * ca / 2 - 1 / 2 * sct * ca + α / 4 * sct * sa

    ofr = θ / 2 * sct * sp
    ofr_θ = (1 / 2 * sct * sp + θ / 2 * cct * sp * θ′_θ) * s * 2
    ofr_ϕ = θ / 2 * sct * cp
    ofr_α = θ / 2 * cct * sp * θ′_α + θ / 2 * sct * cp / 2

    ofi = -θ / 2 * sct * cp
    ofi_θ = (-1 / 2 * sct * cp - θ / 2 * cct * cp * θ′_θ) * s * 2
    ofi_ϕ = θ / 2 * sct * sp
    ofi_α = -θ / 2 * cct * cp * θ′_α + θ / 2 * sct * sp / 2

    return (dgr, dgr_θ, dgr_α,
            dgi, dgi_θ, dgi_α,
            ofr, ofr_θ, ofr_ϕ, ofr_α,
            ofi, ofi_θ, ofi_ϕ, ofi_α)
end
function compute(dri::DetDrive, grad)
    θ = dri.θ
    ϕ = dri.ϕ
    α = dri.α

    (dg1r, dg1r_θ, dg1r_α,
     dg1i, dg1i_θ, dg1i_α,
     of1r, of1r_θ, of1r_ϕ, of1r_α,
     of1i, of1i_θ, of1i_ϕ, of1i_α) = detdrive_terms(θ, ϕ, α, 1)
    (dg3r, dg3r_θ, dg3r_α,
     dg3i, dg3i_θ, dg3i_α,
     of3r, of3r_θ, of3r_ϕ, of3r_α,
     of3i, of3i_θ, of3i_ϕ, of3i_α) = detdrive_terms(θ, ϕ, α, √2)

    grad[1] = @SMatrix[complex(dg1r_θ, dg1i_θ)   complex(of1r_θ, of1i_θ)   0  0
                       complex(-of1r_θ, of1i_θ)  complex(dg1r_θ, -dg1i_θ)  0  0
                       0  0  complex(dg3r_θ, dg3i_θ)   complex(of3r_θ, of3i_θ)
                       0  0  complex(-of3r_θ, of3i_θ)  complex(dg3r_θ, -dg3i_θ)]
    grad[2] = @SMatrix[0                         complex(of1r_ϕ, of1i_ϕ)  0  0
                       complex(-of1r_ϕ, of1i_ϕ)  0                        0  0
                       0  0  0                         complex(of3r_ϕ, of3i_ϕ)
                       0  0  complex(-of3r_ϕ, of3i_ϕ)  0]
    grad[3] = @SMatrix[complex(dg1r_α, dg1i_α)   complex(of1r_α, of1i_α)   0  0
                       complex(-of1r_α, of1i_α)  complex(dg1r_α, -dg1i_α)  0  0
                       0  0  complex(dg3r_α, dg3i_α)   complex(of3r_α, of3i_α)
                       0  0  complex(-of3r_α, of3i_α)  complex(dg3r_α, -dg3i_α)]
    return @SMatrix[complex(dg1r, dg1i)   complex(of1r, of1i)   0  0
                    complex(-of1r, of1i)  complex(dg1r, -dg1i)  0  0
                    0  0  complex(dg3r, dg3i)   complex(of3r, of3i)
                    0  0  complex(-of3r, of3i)  complex(dg3r, -dg3i)]
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

const pm_use_ramp_drive = true

mutable struct PMPulseSeq{N,S,P,OB,RB}
    const s::S # Sequence
    const params::P # Input parameters
    res::Float64 # Output result
    const grads::P # Output gradients

    const op_buff::OB
    const res_buff::RB

    function PMPulseSeq{N}() where N
        ops = ntuple(Val(N + 1)) do i
            if pm_use_ramp_drive
                return i == 1 ? PhaseZ() : DetDrive()
            else
                return i == 1 ? PhaseZ() : Drive()
            end
        end
        s = Sequence{OPType}(ops)
        params = MVector{N + 2,Float64}(undef)
        grads = MVector{N + 2,Float64}(undef)
        if pm_use_ramp_drive
            op_buff = Vector{OPType}(undef,3N + 1)
            res_buff = MVector{3N + 1,Float64}(undef)
        else
            op_buff = Vector{OPType}(undef,2N + 1)
            res_buff = MVector{2N + 1,Float64}(undef)
        end
        return new{N,typeof(s),typeof(params),typeof(op_buff),typeof(res_buff)}(
            s, params, NaN, grads, op_buff, res_buff)
    end
end

@inline function _pm_set_params_buff!(res_buff, params, N)
    @inbounds begin
        res_buff[1] = params[1]
        angle = params[2] / N
        for i in 1:N
            if pm_use_ramp_drive
                res_buff[3 * i - 1] = angle
                res_buff[3 * i] = params[i + 2]
                res_buff[3 * i + 1] = 0
            else
                res_buff[2 * i] = angle
                res_buff[2 * i + 1] = params[i + 2]
            end
        end
    end
end

function update_params!(ps::PMPulseSeq{N}, params) where N
    if !isnan(ps.res) && all(ps.params .== params)
        return
    end
    @assert length(params) == N + 2
    res_buff = ps.res_buff

    _pm_set_params_buff!(res_buff, params, N)
    set_params(ps.s, res_buff)
    op = compute(ps.s, ps.op_buff)
    res = convert_res_grads(op, ps.op_buff, res_buff)

    @inbounds begin
        grads = ps.grads
        grads[1] = res_buff[1]
        angle_grad = 0.0
        for i in 1:N
            if pm_use_ramp_drive
                angle_grad += res_buff[3 * i - 1]
                grads[i + 2] = res_buff[3 * i]
            else
                angle_grad += res_buff[2 * i]
                grads[i + 2] = res_buff[2 * i + 1]
            end
        end
        grads[2] = angle_grad / N
        ps.res = res
        ps.params .= params
    end
    return
end

function optimize_pulse!(ps::PMPulseSeq{N}, init_params; opt_angle=false) where N
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
    @variable(m, total_angle >= 0, start=abs(init_params[2]))
    @variable(m, phases[i=1:N], start=init_params[i + 2])
    obj = @NLexpression(m, 1e-10 + fidelity(global_z, total_angle, phases...))
    if opt_angle
        obj = @NLexpression(m, (total_angle + 1) * obj)
    end
    @NLobjective(m, Min, obj)
    JuMP.optimize!(m)
    params = [value(global_z); value(total_angle); value.(phases)]
    res = res_func(params...)
    return value(obj), res, params
end

mutable struct PMRampSeq{N,S,P,OB}
    const s::S # Sequence
    const params::P # Input parameters
    res::Float64 # Output result
    const grads::P # Output gradients

    const op_buff::OB

    function PMRampSeq{N}() where N
        ops = ntuple(Val(N + 1)) do i
            return i == 1 ? PhaseZ() : DetDrive()
        end
        s = Sequence{OPType}(ops)
        params = MVector{3N + 1,Float64}(undef)
        grads = MVector{3N + 1,Float64}(undef)
        op_buff = Vector{OPType}(undef,3N + 1)
        return new{N,typeof(s),typeof(params),typeof(op_buff)}(
            s, params, NaN, grads, op_buff)
    end
end

function update_params!(ps::PMRampSeq{N}, params) where N
    if !isnan(ps.res) && all(ps.params .== params)
        return
    end
    @assert length(params) == 3N + 1
    ps.res = NaN
    ps.params .= params
    set_params(ps.s, ps.params)
    op = compute(ps.s, ps.op_buff)
    ps.res = convert_res_grads(op, ps.op_buff, ps.grads)
    return
end

function optimize_pulse!(ps::PMRampSeq{N}, init_params; opt_angle=false) where N
    m = Model(NLopt.Optimizer)
    set_attribute(m, "algorithm", :LD_SLSQP)
    function res_func(params...)
        @assert length(params) == 3N + 1
        update_params!(ps, params)
        return ps.res
    end
    function grad_func(g, params...)
        @assert length(params) == 3N + 1
        update_params!(ps, params)
        g .= ps.grads
        return
    end
    register(m, :fidelity, 3N + 1, res_func, grad_func, autodiff=false)
    @variable(m, global_z, start=init_params[1])
    @variable(m, total_angle >= 0, start=abs(init_params[2]))
    @variable(m, phases[i=1:N + 1], start=init_params[i + 2])
    fid_args = []
    push!(fid_args, global_z)
    start_phase = phases[1]
    angle = total_angle / N
    for i in 1:N
        end_phase = phases[i + 1]
        push!(fid_args, angle)
        push!(fid_args, start_phase)
        push!(fid_args, end_phase - start_phase)
        start_phase = end_phase
    end
    res_expr = @NLexpression(m, fidelity(fid_args...))
    obj = @NLexpression(m, 1e-10 + res_expr)
    # @show res_expr
    if opt_angle
        obj = @NLexpression(m, (total_angle + 1) * obj)
    end
    @NLobjective(m, Min, obj)
    JuMP.optimize!(m)
    params = [value(global_z); value(total_angle); value.(phases)]
    return value(obj), value(res_expr), params
end

# Params after output (first number is single qubit Z rotation phase, second number is total rotation angle, third number on are the phase of each of the small rotation steps)
# [-4.082623222999641, 4.541397221247234, -0.03245830600737785, -0.8124489994546822, 0.47202474200536804, -0.5500803261909537, -0.6177047316919362, 0.5411791619982704, -0.2557655640738533, -0.38314193680002606, 0.3118226680374636, 0.5029169562968007, -0.19251679080546885, -0.6830286871769535, -0.9731886766315152, 0.44235641455254676, 0.22454690222112847, 0.6042372041827836, 0.5738136841262153, -0.6144887991414903, 0.3956674052845687, -0.4650998175643587, 0.7887892642909274, 0.11058986980661852, 0.3804075261097273, 1.1650585873558115, 0.1872911377975865, 0.9567278306607799, -0.4947561155344319, 0.6894623797962316, 0.09709192320588979, 0.7557011190748159, 0.8259314495280358, 0.8185719680309294, 0.9349236332175126, 1.3108418952354415, 0.8231074271575507, 1.3398388812930984, 1.4380091588537953, 1.3046589025558661, 1.0254758451006798, 0.5774038163683318, 0.8788491987152173, -0.07075782181316924, -0.5112282797272782, 0.6423404924404528, 0.027051359134414, 0.6404011131017019, -0.13714724496857644, -0.6502383866425321, 0.6573836386387235, -0.1668395082070339, 0.3338638845614989, -0.6685430326154859, -0.21722996509706732, 0.06410340375354018, -1.1336402705064414, -1.1331811610610962, -1.1679161618322564, -1.1513218730649466, -0.6240970802227468, -1.1024089561607755, -1.1884290426757993, 0.09776880469789495, -0.3297167756824757, 0.09422865956652068, 0.43128269240808786, 0.27626678345552685, -0.8704769881183879, -1.0788257414112599, 0.06424796464472106, -0.3954332707109115, 0.4754886348364921, -0.2384744955069249, -0.7489516745171744, -0.8813498815632604, -1.164959134556153, 0.577457335488964, -1.2008233920845337, -0.3581232596135072, -0.9058520582238452, -0.7263428645537944, 0.6411906806732239, -1.0437935842604629, -0.7276694585207621, 0.2210376247854085, 0.07694243354081644, 0.44670224404186903, -0.6675923734961656, -0.7336651023078886, -0.30899172164173766, -0.6278027558000905, -0.6015399266994913, 0.19024846578435994, 0.8059652093292294, 0.9677410073384545, 0.26660466658662235, 0.5748930920838033, 0.5499026608622101, -0.23602798130692976, 1.030334485675284, 1.1253630804302641]

# Unitary produced.
#        1.0-1.25328e-14im  9.3721e-8-5.99741e-7im         0.0+0.0im                  0.0+0.0im
# 5.39887e-7+2.77478e-7im   -0.588955+0.808165im           0.0+0.0im                  0.0+0.0im
#        0.0+0.0im                0.0+0.0im               -1.0-4.88601e-14im  -3.55222e-7+2.77726e-7im
#        0.0+0.0im                0.0+0.0im         1.55589e-7-4.23209e-7im      0.306263+0.951947im

# with some optimization to minimize the total angle, parameters:
# [2.1626271844766385, 3.814265885045716, -0.3027625629523783, -0.23020884795929514, -0.34740265654339847, -0.2880887187590676, -0.20762296030020164, -0.1320722163709499, -0.35229949220585777, -0.3678611744321624, -0.050197804054804875, 0.007176099968260903, -0.012637228237826302, -0.04006373508073409, 0.01906289422699667, 0.06029977469349966, 0.1065820329456576, 0.23903154729495554, 0.30091061448192613, 0.33256399091013683, 0.4087811387539078, 0.4733400964260669, 0.2498423112754017, 0.5397847931689576, 0.5749224769707465, 0.5000992899617491, 0.6476891954858022, 0.6265295208825927, 0.7084878413602845, 0.7513614572131168, 0.6358811444539859, 0.7112379914530147, 0.7075248952152451, 0.7078110375955886, 0.7968828941389735, 0.5591681553371312, 0.7601293317101493, 0.7179334956469292, 0.6965926261867795, 0.6875117301039986, 0.6692384227865857, 0.5908718801402236, 0.5334224159015633, 0.4285687084462007, 0.43522193096646866, 0.4181545596805752, 0.39111001567359494, 0.30152878705012137, 0.276690779557065, 0.18666569181131407, 0.13809058779850875, 0.04028350278823808, 0.010007114297413529, -0.08977764421648565, -0.10807773832981889, -0.17602916302484764, -0.2447281961817002, -0.308628188357479, -0.370657443212529, -0.44074730932530815, -0.49104678712214267, -0.5452177568435226, -0.5794616104126119, -0.6135631461782087, -0.5659906165481701, -0.6868797371650719, -0.6975931539374481, -0.7045128420816205, -0.721367198736763, -0.7314454851167795, -0.6911900540196375, -0.6823117741253837, -0.6826225789128925, -0.7295709882907276, -0.6780432926098079, -0.5978365346455061, -0.6662637659250759, -0.45994343861155257, -0.5409672927133446, -0.476263095185496, -0.3495688337603712, -0.3200338164039557, -0.3682179998988105, -0.24565733362436165, -0.24263049575374615, -0.1869191552470115, -0.21786328768807048, -0.16171009423202237, -0.04601004354156222, -0.09675029559463566, -0.017821370776642687, 0.1043300895318691, 0.04313853858575128, 0.0647267676063336, 0.24735253030638107, 0.2169828602032127, 0.32191297059000795, 0.5109170416608902, 0.25655810252627753, 0.21120105074754567, 0.23262211316192158, 0.5439568628972513]

# Unitary produced.
#        1.0+1.19724e-10im  0.000203466+9.58723e-5im         0.0+0.0im                  0.0+0.0im
# 3.39435e-5-0.000222346im    -0.557881+0.829921im           0.0+0.0im                  0.0+0.0im
#        0.0+0.0im                  0.0+0.0im               -1.0-2.35807e-10im  0.000236742-1.54863e-5im
#        0.0+0.0im                  0.0+0.0im         7.50387e-5+0.000225069im     0.377537+0.925995im

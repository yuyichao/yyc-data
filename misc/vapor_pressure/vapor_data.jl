#!/usr/bin/julia

@kwdef struct VaporPressureData
    Tmelt::Float64
    # Solid phase
    As::Float64
    Bs::Float64
    Cs::Float64 = 0
    Ds::Float64 = 0
    # Liquid phase
    Al::Float64
    Bl::Float64
    Cl::Float64 = 0
    Dl::Float64 = 0
end

# C to Torr
function vapor_pressure_kernel(T, A, B, C, D)
    T += 273.15
    return exp10(5.006 + A + B / T + C * log10(T) + D / T^3) * 0.00750062
end

function vapor_pressure(T, data::VaporPressureData)
    if T > data.Tmelt
        return vapor_pressure_kernel(T, data.Al, data.Bl, data.Cl, data.Dl)
    else
        return vapor_pressure_kernel(T, data.As, data.Bs, data.Cs, data.Ds)
    end
end

const NaData = VaporPressureData(Tmelt=370.944 - 273.15, As=5.298, Bs=-5603,
                                 Al=4.704, Bl=-5377)
const LiData = VaporPressureData(Tmelt=453.65 - 273.15, As=5.667, Bs=-8310,
                                 Al=5.055, Bl=-8023)
const YbData = VaporPressureData(Tmelt=1097 - 273.15, As=9.111, Bs=-8111, Cs=-1.0849,
                                 Al=9.111, Bl=-8111, Cl=-1.0849)
const CaData = VaporPressureData(Tmelt=1115 - 273.15, As=10.127, Bs=-9517, Cs=-1.4030,
                                 Al=10.127, Bl=-9517, Cl=-1.4030)
const InData = VaporPressureData(Tmelt=429.7485 - 273.15, As=5.991, Bs=-12548,
                                 Al=5.374, Bl=-12276)
const TlData = VaporPressureData(Tmelt=577 - 273.15, As=5.971, Bs=-9447,
                                 Al=5.259, Bl=-9037)
const AgData = VaporPressureData(Tmelt=1234.93 - 273.15, As=9.127, Bs=-14999, Cs=0.7845,
                                 Al=5.752, Bl=-13827)
const YData = VaporPressureData(Tmelt=1799 - 273.15, As=9.735, Bs=-22306, Cs=-0.8705,
                                Al=5.795, Bl=-20341)

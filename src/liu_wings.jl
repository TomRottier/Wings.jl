## wing aerofoil and planform implemented from Liu et al. 2006 Avian wing geometry and kinetmatics. AIAA Journal

module LiuWings
using Wings
export LiuPlanform, LiuAerofoil, seagull, merganser, teal

# planform
"""
    LiuPlanform <: AbstractPlanform

Planform as defined in Liu et al. 2006 Avian wing geometry and kinetmatics. AIAA Journal

    LiuPlanform(r, ϕ, E, c₀)

Quarter chord defined by two jointed arm model with two parameters:
- `r`: relative length of upper arm
- `ϕ`: angle between upper and lower arm

Chord distribution defined by 2 parameters:
- `E`: vector of 5 coefficients for their correction function
- `c₀`: root chord length normalised by wing length
"""
struct LiuPlanform <: AbstractPlanform
    r::Float64
    ϕ::Float64
    E::NTuple{5,Float64}
    c₀::Float64
end

# aerofoils
"""
    LiuAerofoil <: AbstractAerofoil

Aerofoil as defined in Liu et al. 2006 Avian wing geometry and kinetmatics. AIAA Journal

    LiuAerofoil(S, A, zcmax, ztmax)

Upper and lower surfaces of aerofoil expressed as a camber line ± thickness distribution
- S: vector of 3 coefficients for Birnbaum-Galuert camber line
- A: vector of 4 coefficients for thickness distribution
- zcmax: vector of 2 coefficients for maximum camber distribution along span
- ztmax: vector of 2 coefficients for maximum thickness distribution along span

Note the owl wing in their paper uses different functions for the camber distributions which is not included here.
"""
struct LiuAerofoil <: AbstractAerofoil
    S::NTuple{3,Float64}
    A::NTuple{4,Float64}
    zcmax::NTuple{2,Float64}
    ztmax::NTuple{2,Float64}
end

# birnbaum-glauert mean camber line, normalised by chord
camber(η, S, zcmax) = zcmax * η * (1 - η) * sum(n -> S[n] * (2η - 1)^(n - 1), 1:3)

# maximum camber distribution along span, normalised by chord
max_camber(ξ, zcmax) = zcmax[1] / (1 + zcmax[2] * ξ^1.4)

# thickness distribution, normalised by chord
thickness(η, A, ztmax) = max(0.001, ztmax * sum(n -> A[n] * (η^(n + 1) - √η), 1:4))

# maximum thickness distribution along span, normalised by chord
max_thickness(ξ, ztmax) = ztmax[1] / (1 + ztmax[2] * ξ^1.4)

# aerofoil value at a single point from upper or lower surface
function liu_aerofoil_pt(η, ξ, S, A, zcmax, ztmax; upper=true)
    __camber = camber(η, S, max_camber(ξ, zcmax))
    __thickness = thickness(η, A, max_thickness(ξ, ztmax))
    surface = upper ? 1.0 : -1.0
    return [η, ξ, __camber + surface * __thickness]
end

# aerofoil
function liu_aerofoil(ξ, S, A, zcmax, ztmax; n=100)
    ηs = Wings.chordwise_coordinates(n ÷ 2)
    upper_surface = map(η -> liu_aerofoil_pt(η, ξ, S, A, zcmax, ztmax; upper=true), ηs)
    lower_surface = map(η -> liu_aerofoil_pt(η, ξ, S, A, zcmax, ztmax; upper=false), ηs)

    return [upper_surface; reverse(lower_surface)]
end

# chord distribution
Fok(ξ) = 0.0 ≤ ξ ≤ 0.5 ? 1.0 : 4ξ * (1 - ξ)
Fcorr(ξ, E) = sum(n -> E[n] * (ξ^(n + 2) - ξ^8), 1:5)
liu_chord(ξ, E, c₀) = c₀ * (Fok(ξ) .+ Fcorr(ξ, E))

# quarter chord line
wing_flexion(ξ, r, ϕ) = ξ < r ? 0.0 : ϕ
liu_quarter_chord(ξ, r, ϕ) = ξ < r ? 0.0 : (ξ - r) * tan(ϕ)

# interface functions
Wings.quarter_chord(ξ, p::LiuPlanform) = liu_quarter_chord(ξ, p.r, p.ϕ)
Wings.chord(ξ, p::LiuPlanform) = liu_chord(ξ, p.E, p.c₀)
Wings.aerofoil(ξ, p::LiuAerofoil; nchord=100) = liu_aerofoil(ξ, p.S, p.A, p.zcmax, p.ztmax; n=nchord)
Wings.aerofoil_pt(η, ξ, p::LiuAerofoil; upper, kwargs...) = liu_aerofoil_pt(η, ξ, p.S, p.A, p.zcmax, p.ztmax; upper, kwargs...)




## data from paper
# seagull
const seagull = (
    S=(3.8735, -0.807, 0.771),
    A=(-15.246, 26.482, -18.975, 4.6232),
    E=(26.08, -209.92, 637.21, -945.068, 695.03),
    zcmax=(0.14, 1.333),
    ztmax=(0.1, 3.546),
    c₀=0.388
)

# merganser
const merganser = (
    S=(3.9385, 0.7466, 1.840),
    A=(-23.1743, 58.3057, -64.3674, 25.7629),
    E=(39.1, -323.8, 978.7, -1417.0, 1001.0),
    zcmax=(0.14, 1.333),
    ztmax=(0.05, 4.0),
    c₀=0.423
)

# teal
const teal = (
    S=(3.9917, -0.3677, 0.0239),
    A=(1.7804, -13.6875, 18.276, -8.279),
    E=(-66.1, 435.6, -1203.0, 1664.1, -1130.2),
    zcmax=(0.11, 4.0),
    ztmax=(0.05, 4.0),
    c₀=0.545
)

# owl - ignore owl uses a different function for zc
const owl = (
    S=(3.9733, -0.8497, -2.723),
    A=(-47.683, 124.5329, -127.0874, 45.876),
    E=(6.3421, -7.5178, -70.9649, 188.0651, -160.1678),
    zcmax=(0.04, 1.8, -0.5), # 0.04*(1 + tanh(1.8ξ - 0.5))
    ztmax=(0.04, 1.78),
    c₀=0.677
)

end
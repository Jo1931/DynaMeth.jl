"""
    struct KineticParameterSeidel

Kinetic parameter for Seidel. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterSeidel <: AbstractParameterStruct
    beta = [-5.001, 26.455, -3.145, 1.5308, -4.4526, 15.615, 1.1064, 0.0, 0.0, 0.0, 0.14969, 0.0, 0.062881, 0.0]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterHybrid

Kinetic parameter for Hybrid. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterHybrid <: AbstractParameterStruct
    beta = [-5.6792660298503055, exp(3.6423241339755763), -3.8637010263367224, exp(2.6317919467787187), -1.7838199223758626, exp(4.084616817214869), exp(-0.5407993134813802), exp(-36.04365338911715), exp(-36.04365338911715), exp(-36.04365338911715), exp(-0.6088952141459639), exp(-36.04365338911715), exp(-2.4280536899946137), exp(-36.04365338911715)]
    gibbs = [0.26464788201008044, 21.379592055290235]
    k = [exp(-8.248523890223751), exp(-4.191229722957382)]
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
    neural = load("data/p_neural.jld2")
end
"""
    struct KineticParameterAutoCatSeidel

Kinetic parameter for Seidel with autoCat step. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterAutoCatSeidel <: AbstractParameterStruct
    beta = [-5.001, 26.455, -3.145, 1.5308, -4.4526, 15.615, 1.1064, 0.0, 0.0, 0.0, 0.14969, 0.0, 0.062881, 0.0, -10, 5.0]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterAutoCat

Kinetic parameter for AutoCat Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterAutoCat <: AbstractParameterStruct
    #beta = [-5.001, 26.455, -3.145, 1.5308, -4.4526, 15.615, -10, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    beta = load_object("data/AutoCat_SS_withoutPhi_GGWGraaf_e8.jld2")[1:end-4]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterAutoCatLumped

Kinetic parameter for Lumped AutoCat Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterAutoCatLumped <: AbstractParameterStruct
    beta = [-7.79606, 26.3361, -3.52772, 23.7252, -3.88219, 10.7081, -10, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterNestler

Kinetic parameter for Nestler Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterNestler <: AbstractParameterStruct
    #beta = [5.441 * 1e-4, -45458, 24.701, -54970, 33.321 * 1e-18, 109959, 8.262 * 1e-6, 6.430 * 1e-14, 119570]
    beta = [-17.9733, 10.4514, -9.4315, 12.6383, -12.6593, -25.2810, 8.262 * 1e-6, -2.8845, -27.4907]

    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterGraaf

Kinetic parameter for Graaf Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterGraaf <: AbstractParameterStruct
    beta = [4.89 * 1e7, -113000, 1.09 * 1e5, -87500, 9.64 * 1e11, -152900, 2.16 * 1e-5, 46800, 7.05 * 1e-7, 61700, 6.37 * 1e-9, 84000]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterGraafRefit

Kinetic parameter for Graaf Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterGraafRefit <: AbstractParameterStruct
    beta = [2.240 * 1e7, -106729, 9.205 * 1e1, -45889, 4.241 * 1e13, -149586, 8.206 * 1e-9, 76594, 1.540 * 1e-3, 14936, 3.818 * 1e-9, 97350]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterBussche

Kinetic parameter for Bussche Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterBussche <: AbstractParameterStruct
    #beta = [1.07 * 1e-10, 36696, 1.22 * 1e5, -94765, 3.45338, 1.578 * 1e-3, 17197, 6.621e-16, 124119]
    beta = [-14.5213, -8.4369, -10.0759, 21.7877, 3.45338, -2.4978, -3.9538, -6.4147, -28.5366]
    #beta = [8.5007, -8.4369, 1.4382, 21.7877, 3453.38, 3.2581, -3.9538, 5.0969, -28.5366]

    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterSlotboom

Kinetic parameter for Slotboom Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterSlotboom <: AbstractParameterStruct
    beta = [7.414 * 1e14, -166000, 1.111 * 1e19, -203700, 126.4, 1.099]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterCampos

Kinetic parameter for Lumped Campos Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterCampos <: AbstractParameterStruct
    beta = [7.414 * 1e14, -166000, 1.111 * 1e19, -203700, 126.4, 1.099]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterBrilman

Kinetic parameter for Lumped Brilman Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterBrilman <: AbstractParameterStruct
    #beta = [4.779 * 1e10, 1.142 * 1e5, 3.529 * 1e14, 1.448 * 1e5, 3.052, 9.695 * 1e2, 3.545 * 1e1]
    beta = [-1.666, 26.2561, 0.205786, 33.2914, 3.052, 969.5, 35.45]

    # beta = [24.5901, 26.2561, 33.4972, 33.2914, 3.052, 969.5, 35.45]
    #beta = [14.8494, 21.1011, 25.6869, 29.5868, 0.130489, 58.5271, 0.0]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterBrilmanFull

Kinetic parameter for Lumped Brilman Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterBrilmanFull <: AbstractParameterStruct
    #beta = [4.779 * 1e10, 1.142 * 1e5, 3.529 * 1e14, 1.448 * 1e5, 3.052, 9.695 * 1e2, 3.545 * 1e1]
    beta = [5.40 * 1e12, 1.22 * 1e5, 3.90 * 1e14, 1.35 * 1e5, 2.97 * 1e-7, 8.95 * 1e-7, 1.23 * 1e2, 1.85 * 1e2, 9.36 * 1e-7, 3.61 * 1e2, 1.36 * 1e2, 4.88 * 1e-7, 9.12, 3.55 * 1e2, 5.61 * 1e-5, 4.44 * 1e-4, 9.55 * 1e-9, 4.5 * 1e-8, 8.13 * 1e-7, 5.03 * 1e2, 6.03 * 1e2, 6.5 * 1e-6, 1.47 * 1e1, 8.56, 7.49 * 1e-7, 6.82 * 1e-7, 8.32, 1.01 * 1e2, 2.33 * 1e1, 3.79 * 1e2]

    # beta = [24.5901, 26.2561, 33.4972, 33.2914, 3.052, 969.5, 35.45]
    #beta = [14.8494, 21.1011, 25.6869, 29.5868, 0.130489, 58.5271, 0.0]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 523.15
    R = 8.314
    phimax = 0.9
end
"""
    struct KineticParameterLinear

Kinetic parameter for Linear Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterLinear <: AbstractParameterStruct
    #beta = [4.779 * 1e10, 1.142 * 1e5, 3.529 * 1e14, 1.448 * 1e5, 3.052, 9.695 * 1e2, 3.545 * 1e1]
    A = zeros(9, 9)

    # beta = [24.5901, 26.2561, 33.4972, 33.2914, 3.052, 969.5, 35.45]
    #beta = [14.8494, 21.1011, 25.6869, 29.5868, 0.130489, 58.5271, 0.0]
    x0 = zeros(9)
    r0 = zeros(3)
end
"""
    struct KineticParameterKoschany

Kinetic parameter for Koschany Kinetic. Fields: beta, gibbs, k, Tmean, R, phimax.
"""
@kwdef mutable struct KineticParameterKoschany <: AbstractParameterStruct
    beta = [3.46e-4, 77.5, 0.5, 22.4, 0.44, -6.2, 0.88, -10]
    gibbs = [0.335747716622, 2.184147818E+1]
    k = [79.174, 0.188] * 1e-4
    Tmean = 550.00
    R = 8.314
    phimax = 0.9
end
function get_beta(p::Tp, ::Type{T}=Float64) where {T} where {Tp<:AbstractParameterStruct}
    return convert.(T, p.beta)
end
function get_gibbs(p::Tp, ::Type{T}=Float64) where {T} where {Tp<:AbstractParameterStruct}
    return convert.(T, p.gibbs)
end
function get_k(p::Tp, ::Type{T}=Float64) where {T} where {Tp<:AbstractParameterStruct}
    return convert.(T, p.k)
end
function get_kinetic_param(p::Tp, ::Type{T}=Float64) where {T} where {Tp<:AbstractParameterStruct}
    nt = (;
        beta=p.beta,
        k=p.k,
        gibbs=p.gibbs
    )

    return NamedTuple([(k, convert.(T, v)) for (k, v) in zip(keys(nt), values(nt))])
end


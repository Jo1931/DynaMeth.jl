struct CatalystStateDynamic <: AbstractKineticComponent end
"""
    (phi::CatalystStateDynamic)(y, T, k, gibbs, param::Ty)

Dynamic catalyst state as in Seidel et al. 2018
"""
function (phi::CatalystStateDynamic)(y, T, k, gibbs, param::Ty) where {Ty}# {Ty<:Union{KineticParameterSeidel,KineticParameterAutoCat,KineticParameterAutoCatSeidel,KineticParameterNestler,KineticParameterKoschany}}
    @unpack R, phimax = param
    K1 = exp(-gibbs[1] * 1e3 / (R * T))
    K2 = exp(-gibbs[2] * 1e3 / (R * T))

    phi = (k[1] * y[3] + k[2] * y[4]) * (phimax - y[7]) - (1 / K1 * k[1] * y[2] + 1 / K2 * k[2] * y[5]) * y[7]
end
struct CatalystStateSteady <: AbstractKineticComponent end
"""
    (phi::CatalystStateSteady)(y, T, k, gibbs, param::Ty)

steady catalyst state as in Seidel et al. 2018
"""
function (phi::CatalystStateSteady)(y, T, k, gibbs, param::Ty) where {Ty}# {Ty<:Union{KineticParameterSeidel,KineticParameterAutoCat,KineticParameterAutoCatSeidel,KineticParameterNestler}}
    @unpack R, phimax = param
    K1 = exp(-gibbs[1] * 1e3 / (R * T))
    K2 = exp(-gibbs[2] * 1e3 / (R * T))
    for i = 2:5
        if y[i] < 0
            y[i] = 0
        end
    end
    gamma = (1 - sqrt(K1 * K2 * (y[4] * y[3] + 1e-6) / (y[5] * y[2] + 1e-6))) / (1 + sqrt(K1 * K2 * (y[4] * y[3] + 1e-6) / (y[5] * y[2] + 1e-6)))
    phi = 1 / 2 * (1 - gamma)
    return 0.1 * (-y[7] + phi)
end
struct CatalystStateNone <: AbstractKineticComponent end
"""
    (phi::CatalystStateNone)(y, T, k, gibbs, param::Ty)

set catalyst state to 0.5
"""
function (phi::CatalystStateNone)(y, T, k, gibbs, param::Ty) where {Ty}
    phi = -y[7] + 0.5
end

struct CatalystStateDynamicJump end
function (phi::CatalystStateDynamicJump)(model, y, T, param::Ty) where {Ty<:Union{KineticParameterSeidel,KineticParameterHybrid}}
    @unpack R, k, gibbs, phimax = param

    K1 = @expression(model, exp(-gibbs[1] * 1e3 / (R * T)))
    K2 = @expression(model, exp(-gibbs[2] * 1e3 / (R * T)))

    phi = @expression(model, (k[1] * y[3] + k[2] * y[4]) * (0.9 - y[7]) - (1 / K1 * k[1] * y[2] + 1 / K2 * k[2] * y[5]) * y[7])
end
struct CatalystSteadyStateJump end
function (phi::CatalystSteadyStateJump)(model, y, T, param)
    phi = @expression(model, 0.0)
end
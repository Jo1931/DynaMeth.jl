struct JacobianLumpedSeidel <: AbstractKineticComponent end
"""
    JacobianLumpedSeidel(y, P, para::KineticParameterSeidel)

Calculate Jacobi-matrix of the surface coverages for Seidel kinetic
"""
function (J::JacobianLumpedSeidel)(y, P, para::Ty) where {Ty<:Union{KineticParameterSeidel,KineticParameterHybrid}}
    @unpack beta = para
    param = beta
    J_dot = zeros(eltype(y), 6, 6)
    J_circle = zeros(eltype(y), 6, 6)
    J_star = zeros(eltype(y), 6, 6)

    f_CH3OH = P * y[1] * 1.0
    f_CO2 = P * y[2] * 1.0
    f_CO = P * y[3] * 1.0
    f_H2 = P * y[4] * 1.0
    f_H2O = P * y[5] * 1.0

    ## Denominators of langmuir isotherms
    theta_dot = 1 / (1 + param[11] * f_CO + param[12] * f_CH3OH + f_CO2 * param[14])
    theta_circle = 1 / (1 + param[7] .* sqrt(abs.(f_H2)))
    theta_star = 1 / (1 + param[9] * f_H2O + param[8] * f_CH3OH + param[13] * f_CO2 + param[10] * param[9] / param[7]^2 * f_H2O / f_H2)

    ## dThtheta_dot_Meoh/dy_i

    J_dot[1, 1] = param[12] * theta_dot - param[12]^2 * P * y[1] * theta_dot^2
    J_dot[1, 2] = -param[14] * param[12] * P * y[1] * theta_dot^2
    J_dot[1, 3] = -param[11] * param[12] * P * y[1] * theta_dot^2

    ## dThtheta_dot_CO2/dy_i

    J_dot[2, 1] = -param[14] * param[12] * P * y[2] * theta_dot^2
    J_dot[2, 2] = param[14] * theta_dot - param[14]^2 * P * y[2] * theta_dot^2
    J_dot[2, 3] = -param[11] * param[14] * P * y[2] * theta_dot^2

    ## dThtheta_dot_CO/dy_i

    J_dot[3, 1] = -param[11] * param[12] * P * y[3] * theta_dot^2
    J_dot[3, 2] = -param[14] * param[11] * P * y[3] * theta_dot^2  ##greska(Tamara)15/05/2019   instead param (12).....param(11)
    J_dot[3, 3] = param[11] * theta_dot - param[11]^2 * P * y[3] * theta_dot^2

    ## dTheta_circle/dy_i

    J_circle[4, 4] = 0.5 * (param[7] * 0.5 * (P * abs(y[4]))^(-0.5)) * theta_circle - 0.5 * param[7]^2 * theta_circle^2   # factor 1/2 comes from the decompostion of H2 to elementary hydrogen

    ## dTheta_star/dy_i

    J_star[1, 1] = param[8] * theta_star - param[8]^2 * P * y[1] * theta_star^2
    J_star[1, 2] = -param[8] * param[13] * P * y[1] * theta_star^2
    J_star[1, 4] = param[8] * param[9] * param[10] / param[7]^2 * y[5] / (P * y[4]^2) * P * y[2] / P * theta_star^2
    J_star[1, 5] = -param[8] * param[9] * (1 + param[10] / param[7]^2 * 1 / (P * y[4])) * P * y[2] * theta_star^2

    J_star[2, 1] = -param[13] * param[8] * P * y[2] * theta_star .^ 2
    J_star[2, 2] = param[13] * theta_star - param[13]^2 * P * y[2] * theta_star^2
    J_star[2, 4] = param[13] * param[9] * param[10] / param[7]^2 * y[5] / (P * y[4]^2) * P * y[2] * theta_star^2
    J_star[2, 5] = -param[13] * param[9] * (1 + param[10] / param[7]^2 * 1 / (P * y[4])) * P * y[2] * theta_star^2

    J_star[4, 1] = -param[7] * param[8] * sqrt(abs(P * y[4])) * theta_star^2
    J_star[4, 2] = -param[7] * param[13] * sqrt(abs(P * y[4])) * theta_star^2
    J_star[4, 4] = 0.5 * param[7] * 1 / sqrt(abs(P * y[4])) * theta_star + param[9] * param[10] / param[7]^2 * (P * y[5] / (P * y[4])^2) * theta_star^2
    J_star[4, 5] = -param[7] * sqrt(abs(P * y[4])) * param[9] * (1 + param[10] / param[7]^2 * 1 / (P * y[4])) * theta_star^2

    J_star[5, 1] = -param[9] * param[8] * P * y[5] * theta_star^2
    J_star[5, 2] = -param[9] * param[13] * P * y[5] * theta_star^2
    J_star[5, 4] = param[9] * param[9] * param[10] / param[7]^2 * P * y[5] / ((P * y[4])^2) * P * y[5] * theta_star^2
    J_star[5, 5] = param[9] * theta_star - param[9]^2 * (1 + param[10] / param[7]^2 * 1 / (P * y[4])) * P * y[5] * theta_star^2

    J = J_dot + J_circle + J_star
end
struct JacobianLumpedNone <: AbstractKineticComponent end
"""
    JacobianLumpedNone(y, P, para)

returns zeros(6,6)
"""
function (J::JacobianLumpedNone)(y, P, para)
    J = zeros(6, 6)
end

struct JacobianLumpedSeidelJump end
function (J::JacobianLumpedSeidelJump)(model, y, Pe, parameter::Ty) where {Ty<:Union{KineticParameterSeidel,KineticParameterHybrid}}
    f_CH3OH = Pe * y[1] * 1.0
    f_CO2 = Pe * y[2] * 1.0
    f_CO = Pe * y[3] * 1.0
    f_H2 = Pe * y[4] * 1.0
    f_H2O = Pe * y[5] * 1.0
    param = parameter.beta

    theta_dot = @expression(model, 1 / (1 + param[11] * f_CO + param[12] * f_CH3OH + f_CO2 * param[14]))
    theta_circle = @expression(model, 1 / (1 + param[7] * sqrt(f_H2)))
    theta_star = @expression(model, 1 / (1 + param[9] * f_H2O + param[8] * f_CH3OH + param[13] * f_CO2 + param[10] * param[9] / param[7]^2 * f_H2O / f_H2))


    ## dThtheta_dot_Meoh/dy_i

    J_dot11 = @expression(model, param[12] * theta_dot - param[12]^2 * Pe * y[1] * theta_dot^2)
    J_dot12 = @expression(model, -param[14] * param[12] * Pe * y[1] * theta_dot^2)
    J_dot13 = @expression(model, -param[11] * param[12] * Pe * y[1] * theta_dot^2)

    ## dThtheta_dot_CO2/dy_i

    J_dot21 = @expression(model, -param[14] * param[12] * Pe * y[2] * theta_dot^2)
    J_dot22 = @expression(model, param[14] * theta_dot - param[14]^2 * Pe * y[2] * theta_dot^2)
    J_dot23 = @expression(model, -param[11] * param[14] * Pe * y[2] * theta_dot^2)

    ## dThtheta_dot_CO/dy_i

    J_dot31 = @expression(model, -param[11] * param[12] * Pe * y[3] * theta_dot^2)
    J_dot32 = @expression(model, -param[14] * param[11] * Pe * y[3] * theta_dot^2)  ##greska(Tamara)15/05/2019   instead param (12).....param(11)
    J_dot33 = @expression(model, param[11] * theta_dot - param[11]^2 * Pe * y[3] * theta_dot^2)

    J_dot = [J_dot11 J_dot12 J_dot13 0 0 0;
        J_dot21 J_dot22 J_dot23 0 0 0;
        J_dot31 J_dot32 J_dot33 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0]
    ## dTheta_circle/dy_i

    J_circle44 = @expression(model, (param[7] * 0.5 * (Pe * y[4])^(-0.5)) * theta_circle - 0.5 * param[7]^2 * theta_circle^2)   # factor 1/2 comes from the decompostion of H2 to elementary hydrogen

    J_circle = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 J_circle44 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0]
    ## dTheta_star/dy_i

    J_star11 = @expression(model, param[8] * theta_star - param[8]^2 * Pe * y[1] * theta_star^2)
    J_star12 = @expression(model, -param[8] * param[13] * Pe * y[1] * theta_star^2)
    J_star14 = @expression(model, param[8] * param[9] * param[10] / param[7]^2 * y[5] / (Pe * y[4]^2) * Pe * y[2] / Pe * theta_star^2)
    J_star15 = @expression(model, -param[8] * param[9] * (1 + param[10] / param[7]^2 * 1 / (Pe * y[4])) * Pe * y[2] * theta_star^2)

    J_star21 = @expression(model, -param[13] * param[8] * Pe * y[2] * theta_star^2)
    J_star22 = @expression(model, param[13] * theta_star - param[13]^2 * Pe * y[2] * theta_star^2)
    J_star24 = @expression(model, param[13] * param[9] * param[10] / param[7]^2 * y[5] / (Pe * y[4]^2) * Pe * y[2] * theta_star^2)
    J_star25 = @expression(model, -param[13] * param[9] * (1 + param[10] / param[7]^2 * 1 / (Pe * y[4])) * Pe * y[2] * theta_star^2)

    J_star41 = @expression(model, -param[7] * param[8] * sqrt(Pe * y[4]) * theta_star^2)
    J_star42 = @expression(model, -param[7] * param[13] * sqrt(Pe * y[4]) * theta_star^2)
    J_star44 = @expression(model, 0.5 * param[7] * 1 / sqrt(Pe * y[4]) * theta_star + param[9] * param[10] / param[7]^2 * (Pe * y[5] / (Pe * y[4])^2) * theta_star^2)
    J_star45 = @expression(model, -param[7] * sqrt(Pe * y[4]) * param[9] * (1 + param[10] / param[7]^2 * 1 / (Pe * y[4])) * theta_star^2)

    J_star51 = @expression(model, -param[9] * param[8] * Pe * y[5] * theta_star^2)
    J_star52 = @expression(model, -param[9] * param[13] * Pe * y[5] * theta_star^2)
    J_star54 = @expression(model, param[9] * param[9] * param[10] / param[7]^2 * Pe * y[5] / ((Pe * y[4])^2) * Pe * y[5] * theta_star^2)
    J_star55 = @expression(model, param[9] * theta_star - param[9]^2 * (1 + param[10] / param[7]^2 * 1 / (Pe * y[4])) * Pe * y[5] * theta_star^2)

    J_star = [J_star11 J_star12 0 J_star14 J_star15 0;
        J_star21 J_star22 0 J_star24 J_star25 0;
        0 0 0 0 0 0;
        J_star41 J_star42 0 J_star44 J_star45 0;
        J_star51 J_star52 0 J_star54 J_star55 0;
        0 0 0 0 0 0]

    J = @expression(model, [i = 1:6, j = 1:6], J_star[i, j] + J_dot[i, j] + J_circle[i, j])
end

struct JacobianLumpedNoneJump end
function (J::JacobianLumpedNoneJump)(model, y, Pe, parameter)

    J_star = [1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1]


    J = @expression(model, [i = 1:6, j = 1:6], J_star[i, j])
end

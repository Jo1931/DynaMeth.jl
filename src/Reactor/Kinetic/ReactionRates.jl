struct ReactionRatesSeidel <: AbstractKineticComponent end
"""
    (r::ReactionRatesSeidel)(y, T, Pe, beta, constants)

Calculate reactions rates for seidel kinetic
"""
function (r::ReactionRatesSeidel)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean

    ## reaction rate
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3 = GGW_const(T)

    ## Partial pressures

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    ## Surface Coverages
    eta_dot = 1 / (1 + beta[11] * f_CO + beta[12] * f_CH3OH + f_CO2 * beta[14])
    eta_circle = 1 / (1 + beta[7] .* abs(f_H2) .^ (0.5))
    eta_star = 1 / (1 + beta[9] .* f_H2O + beta[8] * f_CH3OH + beta[13] * f_CO2 + beta[9] * beta[10] / beta[7]^2 * f_H2O / f_H2)

    ## Reaction Rates
    r1 = k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)
    r2 = k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)
    r3 = k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot)

    ## Considering impact of catalyst dynamics
    # r = [(1 - y[7]) * r1;
    #     (y[7])^(2) * r2;
    #     y[7] / (1 - y[7]) * r3]
    r = [r1; r2; r3]

end
struct ReactionRatesHybrid <: AbstractKineticComponent end
"""
    (r::ReactionRatesSeidel)(y, T, Pe, beta, constants)

Calculate reactions rates for seidel kinetic
"""
function (r::ReactionRatesHybrid)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    p = constants.neural
    ## reaction rate
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3 = GGW_Graaf(T)

    ## Partial pressures

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    ## Surface Coverages
    eta_dot = 1 / (1 + beta[11] * f_CO + beta[12] * f_CH3OH + f_CO2 * beta[14])
    eta_circle = 1 / (1 + beta[7] .* abs(f_H2) .^ (0.5))
    eta_star = 1 / (1 + beta[9] .* f_H2O + beta[8] * f_CH3OH + beta[13] * f_CO2 + beta[9] * beta[10] / beta[7]^2 * f_H2O / f_H2)

    ## Reaction Rates
    r1 = k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)
    r2 = k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)
    r3 = k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot)

    W1 = p.layer_1.weight
    b1 = p.layer_1.bias
    h1 = myswish_0.(W1 * y[1:7] .+ b1)   # oder tanh/sigmoid je nach act

    # Schicht 2
    W2 = p.layer_2.weight
    b2 = p.layer_2.bias
    h2 = myswish_0.(W2 * h1 .+ b2)

    # Ausgabe
    W3 = p.layer_3.weight
    b3 = p.layer_3.bias
    ceta = mysoftplus_0.(W3 * h2 .+ b3)
    ## Considering impact of catalyst dynamics
    r = [ceta[1] * r1;
        ceta[2] * r2;
        ceta[3] * r3]

end
struct ReactionRatesAutoCatSeidel <: AbstractKineticComponent end
"""
    (r::ReactionRatesAutoCatSeidel)(y, T, Pe, beta, constants)

Calculate reactions rates for seidel kinetic with autocatalytic step
"""
function (r::ReactionRatesAutoCatSeidel)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean

    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[15] - beta[16] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3 = GGW_const(T)

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    eta_dot = 1 / (1 + beta[11] * f_CO + beta[12] * f_CH3OH + f_CO2 * beta[14])
    eta_circle = 1 / (1 + beta[7] .* abs(f_H2) .^ (0.5))
    eta_star = 1 / (1 + beta[9] .* f_H2O + beta[8] * f_CH3OH + beta[13] * f_CO2 + beta[9] * beta[10] / beta[7]^2 * f_H2O / f_H2)

    ## Reaction Rates
    r1 = k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)
    r2 = k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)
    r3 = k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot)
    r4 = k4_reac * Pe^2 * y[2] * y[1] * sqrt(Pe * y[4]) * (1 - y[1] * y[5] / (1e-10 + y[2] * y[4]^3 * Pe^2)) * eta_star^2 * eta_circle


    r = [(1 - y[7]) * r1;
        (y[7]) * r2 + y[7]^2 * r4;
        y[7] / (1 - y[7]) * r3]
end

struct ReactionRatesAutoCat <: AbstractKineticComponent end
"""
    (r::ReactionRatesAutoCat)(y, T, Pe, beta, constants)

Calculate reactions rates for AutoCat Kinetic
"""
function (r::ReactionRatesAutoCat)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)

    K_P3 = 1 / K_P3_
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages
    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + 0.0 * beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + 0 * beta[14] * f_CH3OH / f_H2 + 0 * beta[15] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + 0 * beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + 0 * beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + 0 * beta[23] * f_CH3OH / f_H2 + 0 * beta[24] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[25] * f_CH3OH + 0 * beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + 0 * beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + beta[14] * f_CH3OH / f_H2 + beta[15] * f_CH3OH / _protected_sqrt(f_H2) + beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / _protected_sqrt(f_H2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het
    r2 = k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het
    r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    # r3 = k3_reac * (beta[12] *eta_s1 + beta[26]*eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = k4_reac * _protected_sqrt(f_H2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = [r1;
        (r2 + r4);
        r3]


end
function (r::ReactionRatesAutoCat)(y, T, Pe, beta, constants, i::Int)
    Tmean = constants.Tmean
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)

    K_P3 = 1 / K_P3_
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages
    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + 0.0 * beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + 0 * beta[14] * f_CH3OH / f_H2 + 0 * beta[15] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + 0 * beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + 0 * beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + 0 * beta[23] * f_CH3OH / f_H2 + 0 * beta[24] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[25] * f_CH3OH + 0 * beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + 0 * beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + beta[14] * f_CH3OH / f_H2 + beta[15] * f_CH3OH / _protected_sqrt(f_H2) + beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / _protected_sqrt(f_H2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het
    r2 = k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het
    r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    # r3 = k3_reac * (beta[12] *eta_s1 + beta[26]*eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = k4_reac * _protected_sqrt(f_H2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = [r1;
        (r2);
        r3;
        r4]


end
struct ReactionRatesAutoCatRwgs <: AbstractKineticComponent end
"""
    (r::ReactionRatesAutoCat)(y, T, Pe, beta, constants)

Calculate reactions rates for AutoCat Kinetic with different reaction rates for the 2 RWGS centers
"""
function (r::ReactionRatesAutoCatRwgs)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)

    K_P3 = 1 / K_P3_
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages
    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + 0.0 * beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + 0 * beta[14] * f_CH3OH / f_H2 + 0 * beta[15] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + 0 * beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + 0 * beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + 0 * beta[23] * f_CH3OH / f_H2 + 0 * beta[24] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[25] * f_CH3OH + 0 * beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + 0 * beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + beta[14] * f_CH3OH / f_H2 + beta[15] * f_CH3OH / _protected_sqrt(f_H2) + beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / _protected_sqrt(f_H2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het
    r2 = k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het
    #r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r3 = k3_reac * (beta[12] * eta_s1 + beta[26] * eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = k4_reac * _protected_sqrt(f_H2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = [r1;
        (r2 + r4);
        r3]


end
struct ReactionRatesAutoCatR4 <: AbstractKineticComponent end
"""
    (r::ReactionRatesAutoCatR4)(y, T, Pe, beta, constants)

Calculate reactions rates for AutoCat Kinetic only R4
"""
function (r::ReactionRatesAutoCatR4)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)

    K_P3 = 1 / K_P3_
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages
    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + 0.0 * beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + 0 * beta[14] * f_CH3OH / f_H2 + 0 * beta[15] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + 0 * beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + 0 * beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + 0 * beta[23] * f_CH3OH / f_H2 + 0 * beta[24] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[25] * f_CH3OH + 0 * beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + 0 * beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + beta[14] * f_CH3OH / f_H2 + beta[15] * f_CH3OH / _protected_sqrt(f_H2) + beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / _protected_sqrt(f_H2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het
    r2 = k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het
    #r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = k4_reac * _protected_sqrt(f_H2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = [r1;
        r4;
        r3]


end
struct ReactionRatesAutoCatR2 <: AbstractKineticComponent end
"""
    (r::ReactionRatesAutoCatR2)(y, T, Pe, beta, constants)

Calculate reactions rates for AutoCat Kinetic only R2
"""
function (r::ReactionRatesAutoCatR2)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)

    K_P3 = 1 / K_P3_
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages
    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + 0.0 * beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + 0 * beta[14] * f_CH3OH / f_H2 + 0 * beta[15] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + 0 * beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + 0 * beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + 0 * beta[23] * f_CH3OH / f_H2 + 0 * beta[24] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[25] * f_CH3OH + 0 * beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + 0 * beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + beta[14] * f_CH3OH / f_H2 + beta[15] * f_CH3OH / _protected_sqrt(f_H2) + beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / _protected_sqrt(f_H2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2))

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het
    r2 = k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het
    #r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = k4_reac * _protected_sqrt(f_H2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = [r1;
        r2;
        r3]


end
struct ReactionRatesAutoCatWithoutCO <: AbstractKineticComponent end
"""
    (r::ReactionRatesAutoCat)(y, T, Pe, beta, constants)

Calculate reactions rates for AutoCat Kinetic
"""
function (r::ReactionRatesAutoCatWithoutCO)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    #k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)

    K_P3 = 1 / K_P3_
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages

    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / _protected_sqrt(f_H2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = 0 #k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het
    r2 = k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het
    r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    # r3 = k3_reac * (beta[12] *eta_s1 + beta[26]*eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = k4_reac * _protected_sqrt(f_H2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = [r1;
        (r2 + r4);
        r3]


end
struct ReactionRatesAutoCatWithoutCOR2 <: AbstractKineticComponent end
"""
    (r::ReactionRatesAutoCat)(y, T, Pe, beta, constants)

Calculate reactions rates for AutoCat Kinetic
"""
function (r::ReactionRatesAutoCatWithoutCOR2)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    #k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)

    K_P3 = 1 / K_P3_
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages

    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    eta_s2 = 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / _protected_sqrt(f_H2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2))

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = 0 #k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het
    r2 = k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het
    r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    # r3 = k3_reac * (beta[12] *eta_s1 + beta[26]*eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = k4_reac * _protected_sqrt(f_H2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = [r1;
        (r2);
        r3]


end
struct ReactionRatesAutoCatLumped end
"""
    (r::ReactionRatesAutoCatLumped)(y, T, Pe, beta, constants)

Calculate reactions rates for Lumped AutoCat Kinetic
"""
function (r::ReactionRatesAutoCatLumped)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean

    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    k3_reac = exp(beta[5] - beta[6] * (Tmean / T - 1))
    k4_reac = exp(beta[7] - beta[8] * (Tmean / T - 1))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3 = GGW_const(T)

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages
    eta_het = 1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    eta_s1 = 1 / (1 + beta[12] .* f_CO + beta[13] * f_CH3OH + beta[14] * f_CO2)
    eta_s2 = 1 / (1 + beta[15] * f_CO2 + beta[16] * f_CH3OH + beta[17] * f_CO)


    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[17] * eta_s2)
    r1 = k1_reac * (f_CO * f_H2 .^ 2 - f_CH3OH / (K_P1)) * eta_s1 * eta_het^4
    r2 = k2_reac * (f_CO2 * f_H2^2 * _protected_sqrt(f_H2) - f_CH3OH * f_H2O / (_protected_sqrt(f_H2) * K_P2)) * eta_s2 * eta_het^5
    r3 = k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het^2
    r4 = k4_reac * (f_CO2 * f_H2^2 * _protected_sqrt(f_H2) * f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (_protected_sqrt(f_H2) * K_P2)) * eta_s2^2 * eta_het^5

    # r = [(1 - y[7]) * r1;
    #     (y[7]) .^ 2 * r2 + y[7]^2 * r4;
    #     y[7] / (1 - y[7]) * r3]

    r = [r1;
        (r2 + r4);
        r3]


end
struct ReactionRatesNestler end
"""
    (r::ReactionRatesNestler)(y, T, Pe, beta, constants)

Calculate reactions rates for Nestler Kinetic
"""
function (r::ReactionRatesNestler)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    R = 8.314
    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end

    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))

    K1 = exp(beta[5] - beta[6] * (Tmean / T - 1))
    K2 = beta[7]
    K3 = exp(beta[8] - beta[9] * (Tmean / T - 1))

    # k1_reac = beta[1] * exp(beta[2] / R / T)
    # k2_reac = beta[3] * exp(beta[4] / R / T)

    # K1 = beta[5] * exp(beta[6] / R / T)
    # K2 = beta[7]
    # K3 = beta[8] * exp(beta[9] / R / T)


    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2_, K_P3_ = GGW_Graaf(T)

    #Umrechnung in Pa
    K_P2 = (K_P2_) * 1e-10
    K_P3 = (K_P3_)


    f_CH3OH = Pe * y[1] * 1e5
    f_CO2 = Pe * y[2] * 1e5
    f_CO = Pe * y[3] * 1e5
    f_H2 = Pe * y[4] * 1e5
    f_H2O = Pe * y[5] * 1e5
    #_protected_sqrt(x::T) where {T} = real.(sqrt.(Complex.(x)))
    ## Surface Coverages
    EQ1 = f_CO2 * _protected_sqrt(f_H2) * f_H2 - _protected_sqrt(f_H2) * f_H2 * f_CH3OH * f_H2O / (f_H2^3 * K_P2)
    EQ2 = f_CO2 * f_H2 - f_CO * f_H2O / (K_P3)

    r1 = 0
    r2 = k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = k2_reac * K2 * EQ2 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)

    # r = [(1 - y[7]) * r1;
    #     (y[7])^2 * (r2);
    #     y[7] / (1 - y[7]) * r3]
    r = [r1; r2; r3]
end
struct ReactionRatesGraaf end
"""
    (r::ReactionRatesGraaf)(y, T, Pe, beta, constants)

Calculate reactions rates for Graaf Kinetic
"""
function (r::ReactionRatesGraaf)(y, T, Pe, beta, constants)

    R = 8.314

    k1_reac = beta[1] * exp(beta[2] / R / T)
    k2_reac = beta[3] * exp(beta[4] / R / T)
    k3_reac = beta[5] * exp(beta[6] / R / T)

    K1 = beta[5] * exp(beta[6] / R / T)
    K2 = beta[7] * exp(beta[8] / R / T)
    K3 = beta[9] * exp(beta[10] / R / T)


    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2_, K_P3_ = GGW_Graaf(T)

    #Umrechnung in Pa
    K_P2 = (K_P2_)
    K_P3 = (K_P3_)

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]


    r1 = k1_reac * K1 * (f_CO * _protected_sqrt(f_H2) * f_H2 - f_CH3OH / (_protected_sqrt(f_H2) * K_P1)) / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r2 = k2_reac * K2 * (f_CO2 * _protected_sqrt(f_H2) * f_H2 - f_CH3OH * f_H2O / (_protected_sqrt(f_H2) * f_H2 * K_P2)) / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O) #k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = k3_reac * K2 * (f_CO2 * f_H2 - f_H2O * f_CO / K_P3) / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r = [r1; r2; r3]
end
struct ReactionRatesGraafRefit end
"""
    (r::ReactionRatesGraaf)(y, T, Pe, beta, constants)

Calculate reactions rates for Graaf Kinetic
"""
function (r::ReactionRatesGraafRefit)(y, T, Pe, beta, constants)

    R = 8.314

    k1_reac = beta[1] * exp(beta[2] / R / T)
    k2_reac = beta[3] * exp(beta[4] / R / T)
    k3_reac = beta[5] * exp(beta[6] / R / T)

    K1 = beta[5] * exp(beta[6] / R / T)
    K2 = beta[7] * exp(beta[8] / R / T)
    K3 = beta[9] * exp(beta[10] / R / T)


    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2_, K_P3_ = GGW_const(T)

    #Umrechnung in Pa
    K_P2 = (K_P2_)
    K_P3 = (K_P3_)

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    Den = (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r1 = k1_reac * K1 * (f_CO * _protected_sqrt(f_H2) * f_H2 - f_CH3OH / (_protected_sqrt(f_H2) * f_H2 * K_P1)) / Den
    r2 = k2_reac * K2 * (f_CO2 * _protected_sqrt(f_H2) * f_H2 - f_CH3OH * f_H2O / (_protected_sqrt(f_H2) * K_P2)) / Den
    r3 = k3_reac * K2 * (f_CO2 * f_H2 - f_H2O * f_CO / K_P3) / Den
    r = [r1; r2; r3]
end
struct ReactionRatesBussche end
"""
    (r::ReactionRatesBussche)(y, T, Pe, beta, constants)

Calculate reactions rates for Bussche Kinetic
"""
function (r::ReactionRatesBussche)(y, T, Pe, beta, constants)

    R = 8.314

    # k1_reac = beta[1] * exp(beta[2] / R / T)
    # k2_reac = beta[3] * exp(beta[4] / R / T)

    # K1 = beta[5]
    # K2 = beta[6] * exp(beta[7] / R / T)
    # K3 = beta[8] * exp(beta[9] / R / T)
    Tmean = 523.15
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))

    K1 = beta[5]
    K2 = exp(beta[6] - beta[7] * (Tmean / T - 1))
    K3 = exp(beta[8] - beta[9] * (Tmean / T - 1))
    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2_, K_P3_ = GGW_Graaf(T)

    #Umrechnung in Pa
    K_P2 = (K_P2_) * 1e-10
    K_P3 = (K_P3_)

    f_CH3OH = Pe * y[1] * 1e5
    f_CO2 = Pe * y[2] * 1e5
    f_CO = Pe * y[3] * 1e5
    f_H2 = Pe * y[4] * 1e5
    f_H2O = Pe * y[5] * 1e5


    r1 = 0
    r2 = k1_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) / ((1 + K1 * f_H2O / f_H2 + K2 * _protected_sqrt(f_H2) + K3 * f_H2O)^3) #k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = k2_reac * (f_CO2 - f_H2O * f_CO / (f_H2 * K_P3)) / (1 + K1 * f_H2O / f_H2 + K2 * _protected_sqrt(f_H2) + K3 * f_H2O)
    r = [r1; r2; r3]
end
struct ReactionRatesSlotboom end
"""
    (r::ReactionRatesSlotboom)(y, T, Pe, beta, constants)

Calculate reactions rates for Slotboom Kinetic
"""
function (r::ReactionRatesSlotboom)(y, T, Pe, beta, constants)

    R = 8.314

    k1_reac = beta[1] * exp(beta[2] / R / T)
    k2_reac = beta[3] * exp(beta[4] / R / T)
    #k3_reac = beta[5] * exp(beta[6] / R / T)

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3 = GGW_Graaf(T)

    #Umrechnung in Pa

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    theta = (_protected_sqrt(f_H2) * beta[6] + f_H2O * beta[5] + f_CH3OH)^(-1)
    r1 = 0# k1_reac * K1 * (f_CO * _protected_sqrt(f_H2) * f_H2 - f_CH3OH / (_protected_sqrt(f_H2) * K_P1)) / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r2 = k1_reac * (f_CO2 * f_H2^2 - f_CH3OH * f_H2O / (f_H2 * K_P2)) * theta^2 #k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = k2_reac * (f_CO2 * _protected_sqrt(f_H2) - f_H2O * f_CO / (K_P3 * _protected_sqrt(f_H2))) * theta
    r = [r1; r2; r3]
end
struct ReactionRatesCampos end
"""
    (r::ReactionRatesCampos)(y, T, Pe, beta, constants)

Calculate reactions rates for Campos Kinetic
"""
function (r::ReactionRatesCampos)(y, T, Pe, beta, constants)

    R = 8.314


    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_const(T)

    K_P3 = 1 / K_P3_
    #Umrechnung in Pa

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]
    sqrt_h2 = _protected_sqrt(f_H2)

    k1_reac = T^(0.103) * exp(-4632.9 / T - 2.934) * 1.559
    k2_reac = T^(-0.234) * exp(191.2 / T - 7.122) * 1.559
    k3_reac = T^(1.875) * exp(-3561.3 / T - 2.776) * 1.559

    Ka = T^(-0.258) * exp(-17.233 - 6289 / T)
    Kb = T^(-0.498) * exp(-15.637 + 6204.9 / T)
    Kc = exp(-19.031 + 7020.3 / T)
    Kd = T^(-0.756) * exp(-20.480 + 11535.3 / T)
    Ke = T^(-1.234) * exp(-17.288 + 13049.6 / T)
    Kf = T^(0.736) * exp(-33.533 + 9702.4 / T)
    Kg = exp(-7.274 + 1409.6 / T)
    Kh = T^(1.036) * exp(-18.450 + 5390.6 / T)

    t_a = (Kc * f_CO + Kd * sqrt_h2 * f_CO2 + 1)^(-1)
    t_b = (Ke * sqrt_h2 * f_CO2 + Kf * sqrt_h2^(-1) * f_CH3OH + 1)^(-1)
    t_c = (Kg * sqrt_h2 + Kh * sqrt_h2^(-1) * f_H2O + 1)^(-1)

    r1 = k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * t_a * t_c
    r2 = k2_reac * (f_CO2 * _protected_sqrt(f_H2) * f_H2 - f_CH3OH * f_H2O / (_protected_sqrt(f_H2) * f_H2 * K_P2)) * t_b * t_c #k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = k3_reac * (f_CO2 * f_H2O - f_H2 * f_CO2^2 / (K_P3 * f_CO + 1e-10)) * t_c * (t_a * Ka + t_b * Kb)
    r = [r1; r2; r3]
end
struct ReactionRatesBrilman end
"""
    (r::ReactionRatesBrilman)(y, T, Pe, beta, constants)

Calculate reactions rates for Slotboom Kinetic
"""
function (r::ReactionRatesBrilman)(y, T, Pe, beta, constants)
    for i = 1:6
        if y[i] < 0
            y[i] = 0
        end
    end
    R = 8.314

    #k1_reac = beta[1] * exp(-beta[2] / R / T)
    #k2_reac = beta[3] * exp(-beta[4] / R / T)
    Tmean = 523.15
    k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    #k3_reac = beta[5] * exp(beta[6] / R / T)

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)
    K_P3 = K_P3_
    #Umrechnung in Pa

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    theta1 = (1 + beta[5] * f_CO2 * f_H2O + beta[6] * f_CO2 * f_H2O / (f_CO .+ 1e-12) / f_H2)^(-1)
    theta2 = (_protected_sqrt(f_H2) + beta[7] * f_CO)^(-1)
    r1 = 0# k1_reac * K1 * (f_CO * _protected_sqrt(f_H2) * f_H2 - f_CH3OH / (_protected_sqrt(f_H2) * K_P1)) / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r2 = k1_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * theta1 * theta2 #k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = k2_reac * (f_CO2 * _protected_sqrt(f_H2) - f_H2O * f_CO / (K_P3 * _protected_sqrt(f_H2))) * theta1 * theta2
    r = [r1; r2; r3]
end
struct ReactionRatesBrilmanFull end
"""
    (r::ReactionRatesBrilman)(y, T, Pe, beta, constants)

Calculate reactions rates for Slotboom Kinetic
"""
function (r::ReactionRatesBrilmanFull)(y, T, Pe, beta, constants)
    for i = 1:6
        if y[i] < 0
            y[i] = 0
        end
    end
    R = 8.314

    k1_reac = beta[1] * exp(-beta[2] / R / T)
    k2_reac = beta[3] * exp(-beta[4] / R / T)
    #Tmean = 523.15
    #k1_reac = exp(beta[1] - beta[2] * (Tmean / T - 1))
    #k2_reac = exp(beta[3] - beta[4] * (Tmean / T - 1))
    #k3_reac = beta[5] * exp(beta[6] / R / T)

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf(T)
    K_P3 = K_P3_
    #Umrechnung in Pa

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    theta_s = (1 + beta[5] * _protected_sqrt(f_H2) * f_CO2 + beta[6] * f_CH3OH / _protected_sqrt(f_H2) + beta[7] * f_CO2 + beta[8] * f_CO + beta[9] * f_CH3OH + beta[10] * f_H2O + beta[11] * f_H2O / _protected_sqrt(f_H2) + beta[12] * _protected_sqrt(f_H2) * f_CO2 / (f_CO + 1e-5) + beta[13] * f_CO2 / (f_CO + 1e-5) + beta[14] * f_H2O / f_H2)^(-1)
    theta_x = (1 + beta[15] * _protected_sqrt(f_H2) + beta[16] * f_CO2 + beta[17] * f_CO + beta[18] * f_CH3OH + beta[19] * f_H2O + beta[20] * f_H2O / _protected_sqrt(f_H2) + beta[21] * f_H2O / f_H2 + beta[22] * f_CO2 / (f_CO + 1e-5))^(-1)
    theta_p = (1 + beta[23] * f_CO2 * _protected_sqrt(f_H2) + beta[24] * f_CO2 + beta[25] * f_CO + beta[26] * f_CH3OH + beta[27] * f_H2O + beta[28] * f_H2O / _protected_sqrt(f_H2) + beta[30] * f_H2O / f_H2 + beta[29] * f_CO2 / (f_CO + 1e-5))^(-1)
    r1 = 0# k1_reac * K1 * (f_CO * _protected_sqrt(f_H2) * f_H2 - f_CH3OH / (_protected_sqrt(f_H2) * K_P1)) / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r2 = k1_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * theta_s * theta_x #k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = k2_reac * (f_CO2 * _protected_sqrt(f_H2) - f_H2O * f_CO / (K_P3 * _protected_sqrt(f_H2))) * theta_p * theta_x
    r = [r1; r2; r3]
end
struct ReactionRatesKoschany end
"""
    (r::ReactionRatesKoschany)(y, T, Pe, beta, constants)

Calculate reactions rates for Koschany Kinetic
"""
function (r::ReactionRatesKoschany)(y, T, Pe, beta, constants)
    Tmean = constants.Tmean
    R = 8.314

    pCH4 = real(Pe * y[1])
    pCO2 = real(Pe * y[2])
    pH2 = real(Pe * y[4])
    pH2O = real(Pe * y[5])


    k_reac = beta[1] * exp(beta[2] * 1e3 / R * (1 / Tmean - 1 / T))

    KOH = beta[3] * exp(beta[4] * 1e3 / R * (1 / Tmean - 1 / T))
    KH2 = beta[5] * exp(beta[6] * 1e3 / R * (1 / Tmean - 1 / T))
    Kmix = beta[7] * exp(beta[8] * 1e3 / R * (1 / Tmean - 1 / T))

    Keq = 137 * real(T)^(-3.998) * exp(158.7 * 1e3 / R / T)

    theta = real(1 / (1 + KOH * pH2O / (real(pH2)^0.5) + KH2 * abs(real(pH2))^0.5 + Kmix * real(pCO2)^0.5))
    r1 = k_reac * (pCO2)^0.5 * (pH2)^0.5 * (1 - (pCH4 * pH2O^2) / (Keq * pCO2 * pH2^4)) * theta^2
    r = [r1; 0; 0]
end


function GGW_const(T)
    A1 = 13.814
    B1 = 3784.7
    #B1 = 3748.7
    C1 = -9.2833
    D1 = 3.1475
    E1 = 4.2613
    A2 = 1581.7
    B2 = 15.0921
    C2 = -8.7639
    D2 = +2.1105e-3
    E2 = 1.9303e-7
    A3 = 1.2777
    B3 = -2167.0
    C3 = 0.5194
    D3 = -1.037e-3
    E3 = 2.331e-7

    K_P1 = (10^(B1 / T + A1 + log10(T) * C1 + D1 * 1e-3 * T - E1 * 1e-7 * T^2))
    K_P2 = (10^(A2 / T + B2 + log10(T) * C2 + D2 * T - E2 * T^2))
    K_P3 = (10^(B3 / T + A3 + log10(T) * C3 + D3 * T + E3 * T^2))

    return [K_P1, K_P2, K_P3]
end

function GGW_const_Jump(model, T)
    A1 = 13.814
    B1 = 3784.7
    #B1 = 3748.7
    C1 = -9.2833
    D1 = 3.1475
    E1 = 4.2613
    A2 = 1581.7
    B2 = 15.0921
    C2 = -8.7639
    D2 = +2.1105e-3
    E2 = 1.9303e-7
    A3 = 1.2777
    B3 = -2167.0
    C3 = 0.5194
    D3 = -1.037e-3
    E3 = 2.331e-7

    K_P1 = @expression(model, (10^(B1 / T + A1 + log10(T) * C1 + D1 * 1e-3 * T - E1 * 1e-7 * T^2)))
    K_P2 = @expression(model, (10^(A2 / T + B2 + log10(T) * C2 + D2 * T - E2 * T^2)))
    K_P3 = @expression(model, (10^(B3 / T + A3 + log10(T) * C3 + D3 * T + E3 * T^2)))

    return [K_P1, K_P2, K_P3]
end
function GGW_Graaf(T)
    R = 8.314

    c1 = -3.94121 * 1e4
    c2 = -5.41516 * 1e1
    c3 = -5.5642 * 1e-2
    c4 = 2.5760 * 1e-5
    c5 = -7.6594 * 1e-9
    c6 = 1.0161 * 1e-12
    c7 = 1.8429 * 1e1

    a1 = 7.44140 * 1e4
    a2 = 1.89260 * 1e2
    a3 = 3.2443 * 1e-2
    a4 = 7.0432 * 1e-6
    a5 = -5.6053 * 1e-9
    a6 = 1.0344 * 1e-12
    a7 = -6.4364 * 1e1

    K_CO = exp((a1 + a2 * T + a3 * T^2 + a4 * T^3 + a5 * T^4 + a6 * T^5 + a7 * T * log(T)) / R / T)
    K_RWGS = exp((c1 + c2 * T + c3 * T^2 + c4 * T^3 + c5 * T^4 + c6 * T^5 + c7 * T * log(T)) / R / T)
    K_CO2 = K_CO * K_RWGS
    return [K_CO, K_CO2, K_RWGS]
end
function GGW_GraafJump(model, T)
    R = 8.314

    c1 = -3.94121 * 1e4
    c2 = -5.41516 * 1e1
    c3 = -5.5642 * 1e-2
    c4 = 2.5760 * 1e-5
    c5 = -7.6594 * 1e-9
    c6 = 1.0161 * 1e-12
    c7 = 1.8429 * 1e1

    a1 = 7.44140 * 1e4
    a2 = 1.89260 * 1e2
    a3 = 3.2443 * 1e-2
    a4 = 7.0432 * 1e-6
    a5 = -5.6053 * 1e-9
    a6 = 1.0344 * 1e-12
    a7 = -6.4364 * 1e1

    K_CO = @expression(model, exp((a1 + a2 * T + a3 * T^2 + a4 * T^3 + a5 * T^4 + a6 * T^5 + a7 * T * log(T)) / R / T))
    K_RWGS = @expression(model, exp((c1 + c2 * T + c3 * T^2 + c4 * T^3 + c5 * T^4 + c6 * T^5 + c7 * T * log(T)) / R / T))
    K_CO2 = @expression(model, K_CO * K_RWGS)
    return [K_CO, K_CO2, K_RWGS]
end
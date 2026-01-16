struct ReactionRatesSeidelJump end
function (r::ReactionRatesSeidelJump)(model, y, T, Pe, param::KineticParameterSeidel)
    @unpack beta, k, R, Tmean = param

    # temperature in K

    A1 = 13.814
    B1 = 3784.7
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


    k1_reac = @expression(model, exp(beta[1] - beta[2] * (Tmean / T - 1)))
    k2_reac = @expression(model, exp(beta[3] - beta[4] * (Tmean / T - 1)))
    k3_reac = @expression(model, exp(beta[5] - beta[6] * (Tmean / T - 1)))

    K_P1 = @expression(model, (10^(B1 / T + A1 + log10(T) * C1 + D1 * 1e-3 * T - E1 * 1e-7 * T^2)))
    K_P2 = @expression(model, (10^(A2 / T + B2 + log10(T) * C2 + D2 * T - E2 * T^2)))
    K_P3 = @expression(model, (10^(B3 / T + A3 + log10(T) * C3 + D3 * T + E3 * T^2)))

    ## Partial pressures
    f_CH3OH = Pe * y[1] * 1.0
    f_CO2 = Pe * y[2] * 1.0
    f_CO = Pe * y[3] * 1.0
    f_H2 = Pe * y[4] * 1.0
    f_H2O = Pe * y[5] * 1.0

    ## Surface Coverages
    eta_dot = @expression(model, 1 / (1 + beta[11] * f_CO + beta[12] * f_CH3OH + f_CO2 * beta[14]))
    eta_circle = @expression(model, 1 / (1 + beta[7] * f_H2^(0.5)))
    eta_star = @expression(model, 1 / (1 + beta[9] * f_H2O + beta[8] * f_CH3OH + beta[13] * f_CO2 + beta[9] * beta[10] / beta[7]^2 * f_H2O / f_H2))

    ## Reaction Rates

    # r = [@NLexpression(model, (1 - y[7]) * k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)),
    #     @NLexpression(model, (y[7])^(2) * k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)),
    #     @NLexpression(model, y[7] / (1 - y[7]) * k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot))]



    r = @NLexpression(model, [i = 1:3],
        [@expression(model, (1 - y[7]) * k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)),
            @expression(model, (y[7])^(2) * k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)),
            @expression(model, y[7] / (1 - y[7]) * k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot))][i])

    # r = @NLexpressions(model, begin
    #     (1 - y[7]) * k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)
    #     (y[7])^(2) * k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)
    #     y[7] / (1 - y[7]) * k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot)
    # end)
    # ## Considering impact of catalyst dynamics
    # #return r

end
struct ReactionRatesHybridJump end
function (r::ReactionRatesHybridJump)(model, y, T, Pe, param::KineticParameterHybrid)
    @unpack beta, k, R, Tmean, neural = param

    # temperature in K

    A1 = 13.814
    B1 = 3784.7
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


    k1_reac = @expression(model, exp(beta[1] - beta[2] * (Tmean / T - 1)))
    k2_reac = @expression(model, exp(beta[3] - beta[4] * (Tmean / T - 1)))
    k3_reac = @expression(model, exp(beta[5] - beta[6] * (Tmean / T - 1)))

    K_P1 = @expression(model, (10^(B1 / T + A1 + log10(T) * C1 + D1 * 1e-3 * T - E1 * 1e-7 * T^2)))
    K_P2 = @expression(model, (10^(A2 / T + B2 + log10(T) * C2 + D2 * T - E2 * T^2)))
    K_P3 = @expression(model, (10^(B3 / T + A3 + log10(T) * C3 + D3 * T + E3 * T^2)))

    ## Partial pressures
    f_CH3OH = Pe * y[1] * 1.0
    f_CO2 = Pe * y[2] * 1.0
    f_CO = Pe * y[3] * 1.0
    f_H2 = Pe * y[4] * 1.0
    f_H2O = Pe * y[5] * 1.0

    ## Surface Coverages
    eta_dot = @expression(model, 1 / (1 + beta[11] * f_CO + beta[12] * f_CH3OH + f_CO2 * beta[14]))
    eta_circle = @expression(model, 1 / (1 + beta[7] * f_H2^(0.5)))
    eta_star = @expression(model, 1 / (1 + beta[9] * f_H2O + beta[8] * f_CH3OH + beta[13] * f_CO2 + beta[9] * beta[10] / beta[7]^2 * f_H2O / f_H2))

    ## Reaction Rates

    # r = [@NLexpression(model, (1 - y[7]) * k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)),
    #     @NLexpression(model, (y[7])^(2) * k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)),
    #     @NLexpression(model, y[7] / (1 - y[7]) * k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot))]
    p = neural
    h1 = @expression(model, [i = 1:5], myswish_0(sum(p.layer_1.weight[i, j] * y[j] for j in 1:7) + p.layer_1.bias[i]))
    #Layer 2
    h2 = @expression(model, [i = 1:5], myswish_0(sum(p.layer_2.weight[i, j] * h1[j] for j in 1:5) + p.layer_2.bias[i]))
    # Output
    ceta = @expression(model, [k = 1:3], mysoftplus_0(sum(p.layer_3.weight[k, j] * h2[j] for j in 1:5) + p.layer_3.bias[k]))

    r = @NLexpression(model, [i = 1:3],
        [@expression(model, ceta[1] * k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)),
            @expression(model, ceta[2] * k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)),
            @expression(model, ceta[3] * k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot))][i])

    # r = @NLexpressions(model, begin
    #     (1 - y[7]) * k1_reac * (Pe^3 * y[3] * y[4]^2 - (y[1] * Pe) / K_P1) * (eta_dot * eta_circle^4)
    #     (y[7])^(2) * k2_reac * (Pe^3 * y[2] * y[4]^2 - (y[1] * y[5] * Pe) / (K_P2 * y[4])) * (eta_star^2 * eta_circle^4)
    #     y[7] / (1 - y[7]) * k3_reac * (Pe * y[2] - (y[3] * y[5] * Pe) / (K_P3 * y[4])) * (eta_star * eta_dot)
    # end)
    # ## Considering impact of catalyst dynamics
    # #return r

end
struct ReactionRatesKoschanyJump end
function (r::ReactionRatesKoschanyJump)(model, y, T, Pe, param)
    @unpack beta, R, Tmean = param


    pCH4 = @expression(model, Pe * y[1])
    pCO2 = @expression(model, Pe * sqrt((y[2])^2))
    pH2 = @expression(model, Pe * sqrt((y[4])^2))
    pH2O = @expression(model, Pe * y[5])


    k_reac = @expression(model, beta[1] * exp(beta[2] * 1e3 / R * (1 / Tmean - 1 / T)))

    KOH = @expression(model, beta[3] * exp(beta[4] * 1e3 / R * (1 / Tmean - 1 / T)))
    KH2 = @expression(model, beta[5] * exp(beta[6] * 1e3 / R * (1 / Tmean - 1 / T)))
    Kmix = @expression(model, beta[7] * exp(beta[8] * 1e3 / R * (1 / Tmean - 1 / T)))

    Keq = @expression(model, 137 * T^(-3.998) * exp(158.7 * 1e3 / R / T))

    theta = @expression(model, (pH2^0.5) / ((pH2^0.5) + KOH * pH2O + KH2 * pH2^0.5 * (pH2^0.5) + Kmix * pCO2^0.5 * (pH2^0.5) + 1e-4))
    r = @NLexpression(model, [i = 1:3],
        [@expression(model, 1e3 * k_reac * (pCO2)^(0.5) * (pH2)^0.5 / (pCO2 + 1e-4) * (pCO2 - (pCH4 * pH2O^2) / (Keq * pH2^4 + 1e-4)) * theta^2),
            0,
            0][i])

end

struct ReactionRatesBusscheJump end
function (r::ReactionRatesBusscheJump)(model, y, T, Pe, param)
    @unpack beta, R, Tmean = param

    R = 8.314

    Tmean = 523.15
    k1_reac = @expression(model, exp(beta[1] - beta[2] * (Tmean / T - 1)))
    k2_reac = @expression(model, exp(beta[3] - beta[4] * (Tmean / T - 1)))

    K1 = beta[5]
    K2 = @expression(model, exp(beta[6] - beta[7] * (Tmean / T - 1)))
    K3 = @expression(model, exp(beta[8] - beta[9] * (Tmean / T - 1)))
    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2_, K_P3_ = GGW_Graaf_Jump(model, T)

    #K_P1, K_P2_, K_P3_ = (1, 1, 1)
    #Umrechnung in Pa
    K_P2 = (K_P2_) * 1e-10
    K_P3 = (K_P3_)

    f_CH3OH = Pe * y[1] * 1e5
    f_CO2 = Pe * y[2] * 1e5
    f_CO = Pe * y[3] * 1e5
    f_H2 = Pe * y[4] * 1e5
    f_H2O = Pe * y[5] * 1e5


    r1 = 0
    r2 = @expression(model, k1_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) / ((1 + K1 * f_H2O / f_H2 + K2 * (f_H2)^(1 / 2) + K3 * f_H2O)^3))#k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = @expression(model, k2_reac * (f_CO2 - f_H2O * f_CO / (f_H2 * K_P3)) / (1 + K1 * f_H2O / f_H2 + K2 * (f_H2)^(1 / 2) + K3 * f_H2O))
    r = @NLexpression(model, [i = 1:3], [r1; r2; r3][i])

end
struct ReactionRatesNestlerJump end
function (r::ReactionRatesNestlerJump)(model, y, T, Pe, param)
    @unpack beta, R, Tmean = param

    k1_reac = @expression(model, exp(beta[1] - beta[2] * (Tmean / T - 1)))
    k2_reac = @expression(model, exp(beta[3] - beta[4] * (Tmean / T - 1)))

    K1 = @expression(model, exp(beta[5] - beta[6] * (Tmean / T - 1)))
    K2 = beta[7]
    K3 = @expression(model, exp(beta[8] - beta[9] * (Tmean / T - 1)))


    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2_, K_P3_ = GGW_Graaf_Jump(model, T)

    #Umrechnung in Pa
    K_P2 = @expression(model, (K_P2_) * 1e-10)
    K_P3 = (K_P3_)


    f_CH3OH = @expression(model, Pe * y[1] * 1e5)
    f_CO2 = @expression(model, Pe * y[2] * 1e5)
    f_CO = @expression(model, Pe * y[3] * 1e5)
    f_H2 = @expression(model, Pe * y[4] * 1e5)
    f_H2O = @expression(model, Pe * y[5] * 1e5)
    #_protected_sqrt(x::T) where {T} = real.(sqrt.(Complex.(x)))
    ## Surface Coverages
    EQ1 = @expression(model, f_CO2 * (f_H2)^(1 / 2) * f_H2 - (f_H2)^(1 / 2) * f_H2 * f_CH3OH * f_H2O / (f_H2^3 * K_P2))
    EQ2 = @expression(model, f_CO2 * f_H2 - f_CO * f_H2O / (K_P3))

    r1 = 0
    r2 = @expression(model, k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / ((f_H2)^(1 / 2) + K3 * f_H2O))
    r3 = @expression(model, k2_reac * K2 * EQ2 / (1 + K1 * f_CO + K2 * f_CO2) / ((f_H2)^(1 / 2) + K3 * f_H2O))

    # r = [(1 - y[7]) * r1;
    #     (y[7])^2 * (r2);
    #     y[7] / (1 - y[7]) * r3]
    r = @NLexpression(model, [i = 1:3], [r1; r2; r3][i])

end
struct ReactionRatesBrilmanJump end
function (r::ReactionRatesBrilmanJump)(model, y, T, Pe, param)
    @unpack beta, R, Tmean = param


    #k1_reac = beta[1] * exp(-beta[2] / R / T)
    #k2_reac = beta[3] * exp(-beta[4] / R / T)
    Tmean = 523.15
    k1_reac = @expression(model, exp(beta[1] - beta[2] * (Tmean / T - 1)))
    k2_reac = @expression(model, exp(beta[3] - beta[4] * (Tmean / T - 1)))
    #k3_reac = beta[5] * exp(beta[6] / R / T)

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf_Jump(model, T)
    K_P3 = K_P3_
    #Umrechnung in Pa

    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    theta1 = @expression(model, 1 / (1 + beta[5] * f_CO2 * f_H2O + beta[6] * f_CO2 * f_H2O / (f_CO .+ 1e-12) / f_H2))
    theta2 = @expression(model, 1 / ((f_H2)^(1 / 2) + beta[7] * f_CO))
    r1 = 0# k1_reac * K1 * (f_CO * _protected_sqrt(f_H2) * f_H2 - f_CH3OH / (_protected_sqrt(f_H2) * K_P1)) / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r2 = @expression(model, k1_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * theta1 * theta2) #k1_reac * K2 * EQ1 / (1 + K1 * f_CO + K2 * f_CO2) / (_protected_sqrt(f_H2) + K3 * f_H2O)
    r3 = @expression(model, k2_reac * (f_CO2 * (f_H2)^(1 / 2) - f_H2O * f_CO / (K_P3 * (f_H2)^(1 / 2))) * theta1 * theta2)
    r = @NLexpression(model, [i = 1:3], [r1; r2; r3][i])
end
struct ReactionRatesAutoCatJump end
function (r::ReactionRatesAutoCatJump)(model, y, T, Pe, param)
    @unpack beta, R, Tmean = param

    # for i = 1:6
    #     if y[i] < 0
    #         y[i] = 0
    #     end
    # end
    k1_reac = @expression(model, exp(beta[1] - beta[2] * (Tmean / T - 1)))
    k2_reac = @expression(model, exp(beta[3] - beta[4] * (Tmean / T - 1)))
    k3_reac = @expression(model, exp(beta[5] - beta[6] * (Tmean / T - 1)))
    k4_reac = @expression(model, exp(beta[7] - beta[8] * (Tmean / T - 1)))

    ## reaction equilibrium constants Kp = 10^(A/T-B)

    K_P1, K_P2, K_P3_ = GGW_Graaf_Jump(model, T)

    K_P3 = @expression(model, 1 / K_P3_)
    f_CH3OH = Pe * y[1]
    f_CO2 = Pe * y[2]
    f_CO = Pe * y[3]
    f_H2 = Pe * y[4]
    f_H2O = Pe * y[5]

    # ## Surface Coverages
    # eta_het = @expression(model,1 / (1 + beta[9] * _protected_sqrt(f_H2) + beta[10] * f_H2O / _protected_sqrt(f_H2) + beta[11] * f_H2O)
    # eta_s1 = @expression(model,1 / (1 + 0.0 * beta[12] .* f_CO + beta[13] * f_CO * _protected_sqrt(f_H2) + 0 * beta[14] * f_CH3OH / f_H2 + 0 * beta[15] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / _protected_sqrt(f_H2))
    # eta_s2 = @expression(model,1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * _protected_sqrt(f_H2) + 0 * beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + 0 * beta[22] * f_CH3OH * f_H2O / f_H2 / _protected_sqrt(f_H2) + 0 * beta[23] * f_CH3OH / f_H2 + 0 * beta[24] * f_CH3OH / _protected_sqrt(f_H2) + 0 * beta[25] * f_CH3OH + 0 * beta[26] * f_CO + beta[27] * f_CO * f_H2O / _protected_sqrt(f_H2) + 0 * beta[28] * f_CH3OH^2 / f_H2 / f_H2)

    eta_het = @expression(model, 1 / (1 + beta[9] * (f_H2)^(1 / 2) + beta[10] * f_H2O / (f_H2)^(1 / 2) + beta[11] * f_H2O))
    eta_s1 = @expression(model, 1 / (1 + beta[12] .* f_CO + beta[13] * f_CO * (f_H2)^(1 / 2) + beta[14] * f_CH3OH / f_H2 + beta[15] * f_CH3OH / (f_H2)^(1 / 2) + beta[16] * f_CH3OH + beta[17] * f_CO2 + beta[18] * f_CO * f_H2O / (f_H2)^(1 / 2)))
    eta_s2 = @expression(model, 1 / (1 + beta[19] * f_CO2 + beta[20] * f_CO2 * (f_H2)^(1 / 2) + beta[21] * f_CH3OH * f_H2O / f_H2 / f_H2 + beta[22] * f_CH3OH * f_H2O / f_H2 / (f_H2)^(1 / 2) + beta[23] * f_CH3OH / f_H2 + beta[24] * f_CH3OH / (f_H2)^(1 / 2) + beta[25] * f_CH3OH + beta[26] * f_CO + beta[27] * f_CO * f_H2O / (f_H2)^(1 / 2) + beta[28] * f_CH3OH^2 / f_H2 / f_H2))

    #eta_het = 1
    #eta_s1 = 1
    #eta_s2 = 1
    #(beta[12] * eta_s1 + beta[26] * eta_s2)
    r1 = @expression(model, k1_reac * (f_CO * f_H2 - f_CH3OH / (f_H2 * K_P1)) * eta_s1 * eta_het)
    r2 = @expression(model, k2_reac * (f_CO2 * f_H2 - f_CH3OH * f_H2O / (f_H2^2 * K_P2)) * eta_s2 * eta_het)
    r3 = @expression(model, k3_reac * (eta_s1 + eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het)
    # r3 = k3_reac * (beta[12] *eta_s1 + beta[26]*eta_s2) * (f_CO * f_H2O^2 / f_H2 - f_CO2 * f_H2O / K_P3) * eta_het
    r4 = @expression(model, k4_reac * (f_H2)^(1 / 2) * (f_CO2 * f_CH3OH - f_CH3OH^2 * f_H2O / (f_H2^3 * K_P2)) * eta_s2^2 * eta_het)

    # r = [(1 - y[7]) * r1;
    #     (y[7]) * r2 + (y[7]) .^ 2 * r4;
    #     r3]

    r = @NLexpression(model, [i = 1:3], [r1;
        @expression(model, (r2 + r4));
        r3][i])


end
struct ReactionRatesLinearJump end
function (r::ReactionRatesLinearJump)(model, y, T, Pe, param)
    @unpack x0, r0, A = param
    dx = @expression(model, [i = 1:9], y[i] - x0[i])
    r = @expression(model, [i = 1:3], r0[i] + sum(A[i, j] * dx[j] for j = 1:9))
end
function GGW_Graaf_Jump(model, T)
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
    return @expression(model, [i = 1:3], [K_CO, K_CO2, K_RWGS][i])
end
function calc_jacobi_rate(x, cstr; h=1e-6)
    model = Model()
    # Ausgangsraten
    r0 = cstr.kinetic.rates(model, x, x[8], cstr.properties.P * 1e-5, cstr.kinetic.parameter)

    m = length(r0)   # Anzahl Reaktionen (bei dir 3)
    n = length(x)    # Länge von x (bei dir vermutlich 8)


    drdx = Matrix{NonlinearExpression}(undef, m, n)

    for j in 1:n
        xh = copy(x)
        xh[j] += h

        r1 = cstr.kinetic.rates(model, xh, xh[8], cstr.properties.P * 1e-5, cstr.kinetic.parameter)

        for i in 1:m
            drdx[i, j] = @NLexpression(model, (r1[i] - r0[i]) / h)
        end
    end


    return value.(r0), x, value.(drdx)
end

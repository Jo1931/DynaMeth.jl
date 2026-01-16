function calcLHSfpo(x, reac::Reactor)
    # unpack parameter
    for i in eachindex(x)
        if x[i] < 0
            x[i] = 0
        end
    end
    @unpack kinetic, properties = reac
    @unpack ndotin, cp_mix, mcat, P, cp_cat, Vgas, qsat, R = properties
    @unpack parameter, jacobianLumped = kinetic
    T = x[8]
    Pe = P / 1e5
    nG = P * Vgas / R / T
    hads = 1e3

    #Jacobian (Theta)
    J = jacobianLumped(x, Pe, parameter)

    #LHS
    M66 = (nG * Matrix(I, 6, 6) + Pe * qsat * mcat * J)
    M69 = hcat(M66, zeros(6, 3))
    M79 = vcat(M69, [0 0 0 0 0 0 1 0 0])
    LHS8 = nG * cp_mix + mcat * cp_cat
    LHS816 = zeros(eltype(x), 1, 6)
    LHS916 = zeros(eltype(x), 1, 6)
    for j = 1:6
        LHS816[1, j] = hads * mcat * qsat * Pe * sum(J[k, j] for k = 1:6)
        LHS916[1, j] = mcat * qsat * Pe * sum(J[k, j] for k = 1:6)
    end
    LHS8v = hcat(LHS816, 0, LHS8, 0)
    M89 = vcat(M79, LHS8v)
    LHS98 = Vgas * P / R / x[8]^2
    LHS9 = hcat(LHS916, 0, LHS98, 0)
    LHS = vcat(M89, LHS9)
    #RHS

    return LHS
end

function calcLHSexp(y, reac::Reactor)
    # unpack parameter
    @unpack properties, kinetic = reac
    @unpack parameter, jacobianLumped = kinetic
    @unpack P, T, Vgas, R, qsat, mcat, yin, ndotin, activity = properties
    Pe = P / 1e5
    nG = P * Vgas / R / T
    #Jacobian (Theta)
    J = jacobianLumped(y, Pe, parameter)

    #LHS
    M66 = (nG * Matrix(I, 6, 6) + Pe * qsat * mcat * (Matrix(I, 6, 6) - y[1:6] * ones(1, 6)) * J)
    M67 = hcat(M66, zeros(6, 1))
    LHS = vcat(M67, [0 0 0 0 0 0 1])
    return LHS
end
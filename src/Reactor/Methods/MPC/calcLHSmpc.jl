function calcLHSmpc(y, reac::Reactor)
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
    M68 = hcat(M66, zeros(6, 2))
    M78 = vcat(M68, [0 0 0 0 0 0 1 0])
    LHS816 = zeros(eltype(y), 1, 6)
    for j = 1:6
        LHS816[1, j] = mcat * qsat * Pe * sum(J[k, j] for k = 1:6)
    end
    M18 = hcat(LHS816, 0, 0)
    LHS = vcat(M78, M18)
    return LHS
end



function calcLHSJumpFpo(model, y, reac::Reactor)
    if reac.kinetic.rates isa ReactionRatesSeidelJump
        calcLHSJumpFpoSeidel(model, y, reac::Reactor)
    else
        calcLHSJumpFpoSteadyPhi(model, y, reac::Reactor)
    end
end

function calcLHSJumpFpoSeidel(model, y, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack ndotin, cp_mix, k, A, qsat, mcat, H, P, cp_cat, Vgas, reflux, Trec, R = properties
    @unpack parameter, jacobianLumped = kinetic

    T = y[8]
    Pe = P / 1e5

    nG = @expression(model, P * Vgas / (R * T))
    J = jacobianLumped(model, y, Pe, parameter)

    eye = Matrix(I, 6, 6)

    nG_eye = @expression(model, [i = 1:6, j = 1:6], eye[i, j] * nG)



    LHS = Matrix{NonlinearExpr}(undef, 9, 9)
    LHS[1:6, 1:6] = @expression(model, [i = 1:6, j = 1:6], nG_eye[i, j] + Pe * mcat * qsat * J[i, j])
    LHS[7, 7] = 1
    LHS[8, 8] = @expression(model, nG * cp_mix + mcat * cp_cat)
    hads = 1e3
    LHS[8, 1:6] = @expression(model, [[1], j = 1:6], hads * mcat * qsat * Pe * sum(J[k, j] for k = 1:6))
    LHS[9, 1:6] = @expression(model, [[1], j = 1:6], mcat * qsat * Pe * sum(J[k, j] for k = 1:6))
    LHS[9, 8] = @expression(model, Vgas * P / R / y[8]^2)

    # M66 = @expression(model, [i = 1:6, j = 1:6], nG_eye[i, j] + Pe * mcat * qsat * J[i, j])
    # M69 = hcat(M66, zeros(6, 3))
    # #M[1:6,1:6] = (nG*Matrix(I,6,6) + Pe* p["mkat"]*p["q_sat"]*( Matrix(I,6,6) - y[1:6]*ones(1,6))*J);   # component balances
    # M79 = vcat(M69, [0 0 0 0 0 0 1 0 0])

    # LHS8 = @expression(model, nG * cp_mix + mcat * cp_cat)
    # ## adsorptionsenthalpie
    # hads = 1e3
    # LHS816 = @expression(model, [[1], j = 1:6], hads * mcat * qsat * Pe * sum(J[k, j] for k = 1:6))

    # LHS8v = hcat(LHS816, 0, LHS8, 0)
    # M89 = vcat(M79, LHS8v)


    # LHS916 = @expression(model, [[1], j = 1:6], mcat * qsat * Pe * sum(J[k, j] for k = 1:6))
    # LHS98 = @expression(model, Vgas * P / R / y[8]^2)
    # LHS9 = hcat(LHS916, 0, LHS98, 0)
    # LHS = vcat(M89, LHS9)
    return LHS
end
function calcLHSJumpFpoSteadyPhi(model, y, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack ndotin, cp_mix, k, A, qsat, mcat, H, P, cp_cat, Vgas, reflux, Trec, R = properties
    @unpack parameter, jacobianLumped = kinetic

    T = y[8]
    Pe = P / 1e5

    nG = @expression(model, P * Vgas / (R * T))
    #J = jacobianLumped(model, y, Pe, parameter)

    eye = Matrix(I, 6, 6)

    nG_eye = @expression(model, [i = 1:6, j = 1:6], eye[i, j] * nG)



    LHS = Matrix{NonlinearExpr}(undef, 9, 9)
    LHS[1:6, 1:6] = @expression(model, [i = 1:6, j = 1:6], nG_eye[i, j])
    LHS[7, 7] = 1
    LHS[8, 8] = @expression(model, nG * cp_mix + mcat * cp_cat)
    hads = 1e3
    LHS[8, 1:6] = @expression(model, [[1], j = 1:6], 0)
    LHS[9, 1:6] = @expression(model, [[1], j = 1:6], 0)
    LHS[9, 8] = @expression(model, Vgas * P / R / y[8]^2)


    return LHS
end
function calcLHSJumpRplug(model, ident, y, reac::Reactor)
    if reac.kinetic.rates isa ReactionRatesSeidelJump
        calcLHSJumpRplugSeidel(model, ident, y, reac::Reactor)
    else
        calcLHSJumpRplugSteadyPhi(model, ident, y, reac::Reactor)
    end
end
function calcLHSJumpRplugSeidel(model, ident, y, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack ndotin, cp_mix, k, A, qsat, mcat, H, P, cp_cat, Vgas, reflux, Trec, R = properties
    @unpack parameter, jacobianLumped = kinetic

    T = y[8]
    Pe = P / 1e5

    nG = @expression(model, P * Vgas / (R * T))
    J = jacobianLumped(model, y, Pe, parameter)

    eye = Matrix(I, 6, 6)

    nG_eye = @expression(model, [i = 1:6, j = 1:6], eye[i, j] * nG)


    LHS = Matrix{NonlinearExpr}(undef, 9, 9)
    LHS[1:6, 1:6] = @expression(model, [i = 1:6, j = 1:6], nG_eye[i, j] + Pe * mcat * qsat * sum(ident[i, k] * J[k, j] for k = 1:6))
    LHS[7, 7] = 1
    LHS[8, 8] = @expression(model, nG * cp_mix + mcat * cp_cat)
    hads = 1e3
    LHS[8, 1:6] = @expression(model, [[1], j = 1:6], hads * mcat * qsat * Pe * sum(J[k, j] for k = 1:6))
    LHS[9, 1:6] = @expression(model, [[1], j = 1:6], mcat * qsat * Pe * sum(J[k, j] for k = 1:6))
    LHS[9, 8] = @expression(model, Vgas * P / R / y[8]^2)

    # M66 = @expression(model, [i = 1:6, j = 1:6], nG_eye[i, j] + Pe * mcat * qsat * J[i, j])
    # M69 = hcat(M66, zeros(6, 3))
    # #M[1:6,1:6] = (nG*Matrix(I,6,6) + Pe* p["mkat"]*p["q_sat"]*( Matrix(I,6,6) - y[1:6]*ones(1,6))*J);   # component balances
    # M79 = vcat(M69, [0 0 0 0 0 0 1 0 0])

    # LHS8 = @expression(model, nG * cp_mix + mcat * cp_cat)
    # ## adsorptionsenthalpie
    # hads = 1e3
    # LHS816 = @expression(model, [[1], j = 1:6], hads * mcat * qsat * Pe * sum(J[k, j] for k = 1:6))

    # LHS8v = hcat(LHS816, 0, LHS8, 0)
    # M89 = vcat(M79, LHS8v)


    # LHS916 = @expression(model, [[1], j = 1:6], mcat * qsat * Pe * sum(J[k, j] for k = 1:6))
    # LHS98 = @expression(model, Vgas * P / R / y[8]^2)
    # LHS9 = hcat(LHS916, 0, LHS98, 0)
    # LHS = vcat(M89, LHS9)
    return LHS
end
function calcLHSJumpRplugSteadyPhi(model, ident, y, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack ndotin, cp_mix, k, A, qsat, mcat, H, P, cp_cat, Vgas, reflux, Trec, R = properties
    @unpack parameter, jacobianLumped = kinetic

    T = y[8]
    Pe = P / 1e5

    nG = @expression(model, P * Vgas / (R * T))
    #J = jacobianLumped(model, y, Pe, parameter)

    eye = Matrix(I, 6, 6)

    nG_eye = @expression(model, [i = 1:6, j = 1:6], eye[i, j] * nG)


    LHS = Matrix{NonlinearExpr}(undef, 9, 9)
    LHS[1:6, 1:6] = @expression(model, [i = 1:6, j = 1:6], nG_eye[i, j])
    LHS[7, 7] = 1
    LHS[8, 8] = @expression(model, nG * cp_mix + mcat * cp_cat)
    hads = 1e3
    LHS[8, 1:6] = @expression(model, [[1], j = 1:6], 0)
    LHS[9, 1:6] = @expression(model, [[1], j = 1:6], 0)
    LHS[9, 8] = @expression(model, Vgas * P / R / y[8]^2)
    return LHS
end

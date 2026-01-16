function calcLHSJumpMpc(model, y, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack ndotin, cp_mix, k, A, qsat, mcat, H, P, cp_cat, Vgas, reflux, Trec, R, T = properties
    @unpack parameter, jacobianLumped = kinetic


    Pe = P / 1e5

    nG = @expression(model, P * Vgas / (R * T))
    J = jacobianLumped(model, y, Pe, parameter)

    eye = Matrix(I, 6, 6)

    nG_eye = @expression(model, [i = 1:6, j = 1:6], eye[i, j] * nG)



    LHS = Matrix{NonlinearExpr}(undef, 8, 8)
    LHS[1:6, 1:6] = @expression(model, [i = 1:6, j = 1:6], nG_eye[i, j] + Pe * mcat * qsat * J[i, j])
    LHS[7, 7] = 1

    LHS[8, 1:6] = @expression(model, [[1], j = 1:6], mcat * qsat * Pe * sum(J[k, j] for k = 1:6))


    # LHS = vcat(M89, LHS9)
    return LHS
end
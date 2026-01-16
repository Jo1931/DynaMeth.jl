function calcRHSJumpMpc(model, y, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack k, A, mcat, H, P, T = properties
    @unpack parameter, stoich, rates, catalyst = kinetic


    #T = y[8]
    Pe = P / 1e5

    r = rates(model, y, T, Pe, parameter)



    sigma = @NLexpression(model, [i = 1:6], sum(stoich[i, j] * r[j] * reac.properties.activity * reac.methods.para[j] for j = 1:3))

    RHS16 = @NLexpression(model, [i = 1:6], mcat * sigma[i])
    RHS7 = catalyst(model, y, T, parameter)
    RHS8 = @NLexpression(model, mcat * sum(sigma[j] for j = 1:6))


    RHS = [RHS16;  # component balances
        RHS7; #catalyst
        RHS8] #flow

    return RHS
end
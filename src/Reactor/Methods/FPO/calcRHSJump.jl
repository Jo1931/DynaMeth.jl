
function calcRHSJumpFpo(model, y, Tc, reac::Reactor)
    if reac.kinetic.rates isa ReactionRatesSeidelJump
        calcRHSJumpFpoSeidel(model, y, Tc, reac::Reactor)
    else
        calcRHSJumpFpoSteadyPhi(model, y, Tc, reac::Reactor)
    end
end
function calcRHSJumpFpoSeidel(model, y, Tc, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack k, A, mcat, H, P, = properties
    @unpack parameter, stoich, rates, catalyst = kinetic

    T = y[8]
    Pe = P / 1e5

    r = rates(model, y, T, Pe, parameter)



    sigma = @NLexpression(model, [i = 1:6], sum(stoich[i, j] * r[j] * reac.properties.activity for j = 1:3))

    RHS16 = @NLexpression(model, [i = 1:6], mcat * sigma[i])
    RHS7 = catalyst(model, y, T, parameter)
    RHS8 = @NLexpression(model, -k * A * (y[8] - Tc) - mcat * sum(H[i] * r[i] * reac.properties.activity for i = 1:3))
    RHS9 = @NLexpression(model, mcat * sum(sigma[j] for j = 1:6))


    RHS = [RHS16;  # component balances
        RHS7; #catalyst
        RHS8;#temperarture
        RHS9] #flow rate];
    return RHS
end
function calcRHSJumpFpoSteadyPhi(model, y, Tc, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack k, A, mcat, H, P, = properties
    @unpack parameter, stoich, rates, catalyst = kinetic

    T = y[8]
    Pe = P / 1e5

    r = rates(model, y, T, Pe, parameter)



    sigma = @NLexpression(model, [i = 1:6], sum(stoich[i, j] * r[j] * reac.properties.activity for j = 1:3))

    RHS16 = @NLexpression(model, [i = 1:6], mcat * sigma[i])
    RHS7 = 0
    RHS8 = @NLexpression(model, -k * A * (y[8] - Tc) - mcat * sum(H[i] * r[i] * reac.properties.activity for i = 1:3))
    RHS9 = @NLexpression(model, mcat * sum(sigma[j] for j = 1:6))


    RHS = [RHS16;  # component balances
        RHS7; #catalyst
        RHS8;#temperarture
        RHS9] #flow rate];
    return RHS
end

function calcRHSJumpRplug(model, ident, y, Tc, H, reac::Reactor)
    if reac.kinetic.rates isa ReactionRatesSeidelJump
        calcRHSJumpRplugSeidel(model, ident, y, Tc, H, reac::Reactor)
    else
        calcRHSJumpRplugSteadyPhi(model, ident, y, Tc, H, reac::Reactor)
    end
end
function calcRHSJumpRplugSeidel(model, ident, y, Tc, H, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack k, A, mcat, P, = properties
    @unpack parameter, stoich, rates, catalyst = kinetic
    T = y[8]
    Pe = P / 1e5

    r = rates(model, y, T, Pe, parameter)
    #ident = @expression(model, [i = 1:6, j = 1:6], eye[i, j] - y[i])

    sigma = @NLexpression(model, [i = 1:6], sum(stoich[i, j] * r[j] * reac.properties.activity for j = 1:3))

    RHS16 = @NLexpression(model, [i = 1:6], mcat * sum(ident[i, j] * sigma[j] for j = 1:6))
    RHS7 = catalyst(model, y, T, parameter)
    RHS8 = @NLexpression(model, -k * A * (T - Tc) - mcat * sum(H[i] * r[i] * reac.properties.activity for i = 1:3))
    RHS9 = @NLexpression(model, mcat * sum(sigma[j] for j = 1:6))

    RHS = [RHS16;  # component balances
        RHS7; #catalyst
        RHS8;#temperarture
        RHS9] #flow rate];
end
function calcRHSJumpRplugSteadyPhi(model, ident, y, Tc, H, reac::Reactor)
    @unpack kinetic, properties, methods = reac
    @unpack k, A, mcat, P, = properties
    @unpack parameter, stoich, rates, catalyst = kinetic
    T = y[8]
    Pe = P / 1e5

    r = rates(model, y, T, Pe, parameter)
    #ident = @expression(model, [i = 1:6, j = 1:6], eye[i, j] - y[i])

    sigma = @NLexpression(model, [i = 1:6], sum(stoich[i, j] * r[j] * reac.properties.activity for j = 1:3))

    RHS16 = @NLexpression(model, [i = 1:6], mcat * sum(ident[i, j] * sigma[j] for j = 1:6))
    #RHS7 = catalyst(model, y, T, parameter)
    RHS8 = @NLexpression(model, -k * A * (T - Tc) - mcat * sum(H[i] * r[i] * reac.properties.activity for i = 1:3))
    RHS9 = @NLexpression(model, mcat * sum(sigma[j] for j = 1:6))

    RHS = [RHS16;  # component balances
        0; #catalyst
        RHS8;#temperarture
        RHS9] #flow rate];
end
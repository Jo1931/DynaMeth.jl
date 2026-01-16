function setRplugCons(model, xt, xz, y, Tc, reac::Reactor)
    H, cp_mix = enthalpie(model, y)
    eye = Matrix(I, 6, 6)
    ident = @expression(model, [i = 1:6, j = 1:6], eye[i, j] - y[i])

    RHS = calcRHSJumpRplug(model, ident, y, Tc, H, reac)
    LHS = calcLHSJumpRplug(model, ident, y, reac)

    if reac.methods.signal(0.1) == 0.0
        ss = 0
    else
        ss = 1
    end

    ss = 1
    @NLconstraint(model, ss * sum(LHS[1, i] * xt[i] for i = 1:6) + reac.properties.len * y[9] * xz[1] == RHS[1])
    @NLconstraint(model, ss * sum(LHS[2, i] * xt[i] for i = 1:6) + reac.properties.len * y[9] * xz[2] == RHS[2])
    @NLconstraint(model, ss * sum(LHS[3, i] * xt[i] for i = 1:6) + reac.properties.len * y[9] * xz[3] == RHS[3])
    @NLconstraint(model, ss * sum(LHS[4, i] * xt[i] for i = 1:6) + reac.properties.len * y[9] * xz[4] == RHS[4])
    @NLconstraint(model, ss * sum(LHS[5, i] * xt[i] for i = 1:6) + reac.properties.len * y[9] * xz[5] == RHS[5])
    if !reac.methods.inert
        @constraint(model, y[6] == 0)
    else
        @NLconstraint(model, ss * sum(LHS[6, i] * xt[i] for i = 1:6) + reac.properties.len * y[9] * xz[6] == RHS[6])
    end
    @NLconstraint(model, ss * xt[7] == RHS[7])
    if reac.methods.isothermal
        @NLconstraint(model, xt[8] + xz[8] == 0)
    else
        @NLconstraint(model, (ss * sum(LHS[8, i] * xt[i] for i = vcat(1:6, 8)) + reac.properties.len * cp_mix * y[9] * xz[8]) / 5000 == RHS[8] / 5000)
    end
    @NLconstraint(model, ss * sum(LHS[9, i] * xt[i] for i = vcat(1:6)) + reac.properties.len * xz[9] == RHS[9])

end

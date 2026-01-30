function setCstrCons(model, xt, y, y_i, reac::Reactor)

    n_0 = y_i[9]
    Tc = y_i[10]
    reflux = reac.properties.reflux
    Trec = reac.properties.Trec

    RHS = calcRHSJumpFpo(model, y, Tc, reac)
    LHS = calcLHSJumpFpo(model, y, reac)

    if reac.methods.signal(0.1) == 0.0
        ss = 0
    else
        ss = 1
    end


    yup = @expression(model, [i = 1:6], [0, 1, 1, 1, 0, 1][i] * y[i])
    nrec = @expression(model, [i = 1:6], y[9] * yup[i] * reflux)


    @NLconstraint(model, ss * sum(LHS[1, i] * xt[i] for i = 1:6) == n_0 * y_i[1] + nrec[1] - y[9] * y[1] + RHS[1])
    @NLconstraint(model, ss * sum(LHS[2, i] * xt[i] for i = 1:6) == n_0 * y_i[2] + nrec[2] - y[9] * y[2] + RHS[2])
    @NLconstraint(model, ss * sum(LHS[3, i] * xt[i] for i = 1:6) == n_0 * y_i[3] + nrec[3] - y[9] * y[3] + RHS[3])
    @NLconstraint(model, ss * sum(LHS[4, i] * xt[i] for i = 1:6) == n_0 * y_i[4] + nrec[4] - y[9] * y[4] + RHS[4])
    @NLconstraint(model, ss * sum(LHS[5, i] * xt[i] for i = 1:6) == n_0 * y_i[5] + nrec[5] - y[9] * y[5] + RHS[5])
    if reac.methods.inert
        @NLconstraint(model, ss * sum(LHS[6, i] * xt[i] for i = 1:6) == n_0 * y_i[6] + nrec[6] - y[9] * y[6] + RHS[6])
    else
        @constraint(model, y[6] == 0)
    end

    @NLconstraint(model, ss * xt[7] == RHS[7])

    if reac.methods.isothermal
        @NLconstraint(model, y[8] == Tc)
    else
        @NLconstraint(model, ss * sum(LHS[8, i] * xt[i] for i = vcat(1:6, 8)) == n_0 * reac.properties.cp_mix * (y_i[8] - y[8]) + sum(nrec[i] for i = 1:6) * reac.properties.cp_mix * (Trec - y[8]) + RHS[8])
    end
    @NLconstraint(model, sum(LHS[9, i] * xt[i] for i = 1:6) + y[9] == n_0 + sum(nrec[i] for i = 1:6) + RHS[9])



end
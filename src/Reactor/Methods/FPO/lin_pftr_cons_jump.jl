function pftrcons1(model, xt, xz, y, n_0, l)
    # A -> C  # 

    k1 = 1.0
    k2 = 0
    r1 = @expression(model, k1 * y[1])
    r2 = @expression(model, k2 * y[3])
    mcat = log(2)
    @NLconstraint(model, xt[1] + n_0 * l * xz[1] == mcat * (-r1 + r2))
    @NLconstraint(model, xt[2] + n_0 * l * xz[2] == 0)
    @NLconstraint(model, xt[3] + n_0 * l * xz[3] == mcat * (r1 - r2))
end





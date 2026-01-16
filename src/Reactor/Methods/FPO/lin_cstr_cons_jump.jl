function cstrcons1(model, xt, y, y_i)
    # A -> C  # 

    k1 = 1.0
    k2 = 0
    r1 = @expression(model, k1 * y[1])
    r2 = @expression(model, k2 * y[3])

    n_0 = y_i[4]
    @NLconstraint(model, xt[1] == n_0 * y_i[1] - n_0 * y[1] - r1 + r2)
    @NLconstraint(model, xt[2] == n_0 * y_i[2] - n_0 * y[2])
    @NLconstraint(model, xt[3] == n_0 * y_i[3] - n_0 * y[3] + r1 - r2)
end
function cstrcons2(model, xt, y, y_i)
    # A -> C  # 
    # B -> C 
    k1 = 1.0
    k2 = 0.1
    r1 = @expression(model, k1 * y[1])
    r2 = @expression(model, k2 * y[2])
    #r1m = @expression(model, 2 * y[3])
    n_0 = y_i[4]
    @NLconstraint(model, xt[1] == n_0 * y_i[1] - n_0 * y[1] - r1)
    @NLconstraint(model, xt[2] == n_0 * y_i[2] - n_0 * y[2] - r2)
    @NLconstraint(model, xt[3] == n_0 * y_i[3] - n_0 * y[3] + r1 + r2)
    # @NLconstraint(model, xt[4] == n_0 * y_i[4] - n_0 * y[4])

end

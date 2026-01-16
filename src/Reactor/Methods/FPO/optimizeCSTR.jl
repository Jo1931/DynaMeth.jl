function optimizeCSTR(init, epsilon, cstr; sca=100.0, max_iter=300)
    # x_ub = setUpperBound(cstr)
    # x_lb = setLowerBound(cstr)
    # u_ub = [1, 1, 1, 1, 1, 1, 1, 1.1, 1, 753.15, 753.15, 5, 1.0, 0.0, 1.0]
    # u_lb = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.005, -1.1, 0.00, 300, 453.15, 0.0, 0.0, 0.0, 0.0]
    # u0 = [0.05316199332175653 0.3052336287810913 0.5121880116125959 0.4801406907619201 * 1 0.41910463843764845 * 1 0.9884748992878516 * 1 0.015967694441066885 * 1 0.7305239171611769 0.12941636628455633 523.15 523.15 4.3 0 0 0]
    # x0 = zeros(9, cstr.fpo.nall, cstr.fpo.mall)
    # x0[:, 1:cstr.fpo.nall, 1] .= vcat(cstr.properties.yin, 0.7, cstr.properties.Tc, cstr.properties.ndotin)


    n = cstr.methods.nall
    m = cstr.methods.mall
    n0 = cstr.methods.n
    ordn = cstr.methods.ordn
    ordm = cstr.methods.ordm

    ti = LinRange(0, 1, n0)
    zi = LinRange(0, 1, 2)

    x0 = init[:x]
    u0 = init[:u]
    x_ub = init[:xub]
    x_lb = init[:xlb]
    u_ub = init[:uub]
    u_lb = init[:ulb]

    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 60.0 * 60, "max_iter" => max_iter,
            "nlp_scaling_obj_target_gradient" => sca,
            "tol" => 1e-12,
            #"bound_relax_factor" => 0.0
        ), add_bridges=false)#,"acceptable_constr_viol_tol" => 1e-18))
    ## define Variables
    set_string_names_on_creation(model, false)

    @variable(model, x_lb[c, i, 1] <= x[c=1:9, i=1:cstr.methods.nall] <= x_ub[c, i, 1], start = x0[c, i, 1])
    @variable(model, u_lb[i] <= u[i=1:15] <= u_ub[i], start = u0[i])

    ## define grid (Gauß Lobatto points)
    tspa, zspa = GL_nodes(ti, zi, ordn, ordm)
    @expression(model, tspan[i=1:n], tspa[i])
    @expression(model, zspan[i=1:m], zspa[i])



    # define period and frequency
    @expression(model, tend, u[7] * 3600)
    @expression(model, zend, 1)
    @expression(model, wd, 2 * pi / tend)

    ## define derivative


    #return res
    dxt = dxdt(model, x, tspan, tend, ordn, n, m)

    ## periodic input function
    per_f = cstr.methods.signal
    register(model, :per_f, 1, per_f; autodiff=true)
    @expression(model, per_shift[i=1:n], per_f(2 * pi * tspan[i] + u[8] * pi))
    @expression(model, per[i=1:n], per_f(2 * pi * tspan[i]))

    #Inlet Expressions
    #xin = @NLexpression(model, [c = 1:10, i = 1:n], 0.0)
    @expression(model, xin1[i=1:n], 0)
    @expression(model, xin2[i=1:n], u[1] * (1 - u[13] * per[i]))
    @expression(model, xin3[i=1:n], u[2] * (1 + u[5] * per[i]))
    @expression(model, xin4[i=1:n], u[3] * (1 - u[14] * per[i]))
    @expression(model, xin5[i=1:n], 0)
    @expression(model, xin6[i=1:n], u[9] * (1 - u[6] * per[i]))
    #@expression(model, xin7[i=1:n], x[7, i]) #Egal für optimierung
    @expression(model, xin8[i=1:n], u[10] * (1 + u[15] * per[i]))
    @expression(model, xin9[i=1:n], u[12] * (1 + u[4] * per_shift[i]))
    @expression(model, xin10[i=1:n], u[11] * (1 + u[15] * per[i]))

    ## CSTR-Equations
    for c = 1:9
        @constraint(model, x[c, 1] == x[c, end])
    end
    for i = 2:cstr.methods.nall
        xin = [xin1[i], xin2[i], xin3[i], xin4[i], xin5[i], xin6[i], 0, xin8[i], xin9[i], xin10[i]]
        setCstrCons(model, dxt[:, i, 1], x[:, i, 1], xin, cstr::Reactor)
    end

    ## Imlet Constraints
    @expression(model, mean_N2_in, sum((xin6[j] * xin9[j] * 0.5 + xin6[j-1] * xin9[j-1] * 0.5) * (tspan[j] - tspan[j-1]) for j = 2:n))
    setInputCons(model, u, mean_N2_in, u0, cstr, cstr.methods.optimize)
    n_in = xin9
    if cstr.methods.setMfcCons
        for i = 1:n
            if cstr.properties.P == 3.0e6
                @NLconstraint(model, n_in[i] * xin6[i] >= 0.000277649)
                @NLconstraint(model, n_in[i] * xin6[i] <= 0.002161695)
            elseif cstr.properties.P == 4.0e6
                @NLconstraint(model, n_in[i] * xin6[i] >= 0.000337145)
                @NLconstraint(model, n_in[i] * xin6[i] <= 0.002221191)
            elseif cstr.properties.P == 5.0e6
                @NLconstraint(model, n_in[i] * xin6[i] >= 0.00039664)
                @NLconstraint(model, n_in[i] * xin6[i] <= 0.00247901)
            elseif cstr.properties.P == 6.0e6
                @NLconstraint(model, n_in[i] * xin6[i] >= 0.000456137)
                @NLconstraint(model, n_in[i] * xin6[i] <= 0.002340183)
            end
            @NLconstraint(model, n_in[i] * xin2[i] >= 6.3177e-5)
            @NLconstraint(model, n_in[i] * xin2[i] <= 0.00126234)

            @NLconstraint(model, n_in[i] * xin3[i] >= 0.00014876)
            @NLconstraint(model, n_in[i] * xin3[i] <= 0.00198342)

            @NLconstraint(model, n_in[i] * xin4[i] >= 0.00022294)
            @NLconstraint(model, n_in[i] * xin4[i] <= 0.0041336)
        end
    end
    # for i = 1:n
    #     @NLconstraint(model, n_in[i] >= cstr.properties.ndotin * 0.1)
    # end
    ## Objective Functions
    n_in = xin9
    y_in2 = xin2
    y_in3 = xin3
    n_in2 = @expression(model, sum((n_in[j] * y_in2[j] + n_in[j-1] * y_in2[j-1]) * 0.5 * (tspan[j] - tspan[j-1]) for j = 2:n))
    n_in3 = @expression(model, sum((n_in[j] * y_in3[j] + n_in[j-1] * y_in3[j-1]) * 0.5 * (tspan[j] - tspan[j-1]) for j = 2:n))

    cstr.methods.obj1(model, x, tspan, n, cstr.properties.mcat)
    cstr.methods.obj2(model, x, n_in2, n_in3, tspan, n, cstr.properties.mcat)
    obj1 = model[:obj1]
    obj2 = model[:obj2]

    if cstr.methods.setObjective == 1
        @NLobjective(model, Max, obj1)
        @NLconstraint(model, obj2 >= epsilon)
    elseif cstr.methods.setObjective == 2
        @NLobjective(model, Max, obj2)
        @NLconstraint(model, obj2 >= epsilon)
    end
    # Solve
    @time solu = JuMP.optimize!(model)

    return model
end
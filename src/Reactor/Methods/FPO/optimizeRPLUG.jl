function optimizeRPLUG!(init, epsilon, iter, rplug; sca=100.0, linear_solver="mumps", max_iter=300)
    # x_ub = setUpperBound(rplug)
    # x_lb = setLowerBound(rplug)
    # u_ub = [1, 1, 1, 1, 1, 1, 1, 1.1, 1, 753.15, 753.15, 5, 1.0, 0.0, 1.0]
    # u_lb = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.005, -1.1, 0.00, 300, 453.15, 0.0, 0.0, 0.0, 0.0]
    # u0 = [0.05316199332175653 0.3052336287810913 0.5121880116125959 0.4801406907619201 * 1 0.41910463843764845 * 1 0.9884748992878516 * 1 0.015967694441066885 * 1 0.7305239171611769 0.12941636628455633 523.15 523.15 4.3 0 0 0]
    # x0 = zeros(9, rplug.fpo.nall, rplug.fpo.mall)
    # x0[:, 1:rplug.fpo.nall, 1] .= vcat(rplug.properties.yin, 0.7, rplug.properties.Tc, rplug.properties.ndotin)


    n = rplug.methods.nall
    m = rplug.methods.mall
    n0 = rplug.methods.n
    m0 = rplug.methods.m
    ordn = rplug.methods.ordn
    ordm = rplug.methods.ordm

    ti = LinRange(0, 1, n0)
    zi = LinRange(0, 1, m0)

    x0 = init[:x]
    u0 = init[:u]
    x_ub = init[:xub]
    x_lb = init[:xlb]
    u_ub = init[:uub]
    u_lb = init[:ulb]
    print(u_lb)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 60.0 * 60 * 1, "max_iter" => max_iter,
            "nlp_scaling_obj_target_gradient" => sca,
            "tol" => 1e-12,
            "linear_solver" => linear_solver,
            #"warm_start_init_point" => "yes",
            #"mu_init" => 1e-6,                      # wichtig für schnellen Warmstart
            #"warm_start_bound_push" => 1e-6,
            #"warm_start_bound_frac" => 1e-6,
            #"warm_start_slack_bound_frac" => 1e-6,
            #"warm_start_slack_bound_push" => 1e-6,
            #"warm_start_mult_bound_push" => 1e-6,
            #"bound_relax_factor" => 0.0
        ), add_bridges=false)#,"acceptable_constr_viol_tol" => 1e-18))
    ## define Variables
    set_string_names_on_creation(model, false)

    @variable(model, x_lb[c, i, j] <= x[c=1:9, i=1:n, j=1:m] <= x_ub[c, i, j], start = x0[c, i, j])
    @variable(model, u_lb[i] <= u[i=1:15] <= u_ub[i], start = u0[i])

    ## define grid (Gauß Lobatto points)
    tspa, zspa = GL_nodes(ti, zi, ordn, ordm)
    @expression(model, tspan[i=1:n], tspa[i])
    @expression(model, zspan[i=1:m], zspa[i])


    # define period and frequency
    @expression(model, tend, u[7] * 3600)
    @expression(model, wd, 2 * pi / tend)
    @expression(model, zend, rplug.properties.len)
    ## define derivative


    #return res
    dxt = dxdt(model, x, tspan, tend, ordn, n, m)
    dxz = dxdz(model, x, zspan, zend, ordm, n, m)
    ## periodic input function
    per_f = rplug.methods.signal
    register(model, :per_f, 1, per_f; autodiff=true)
    @expression(model, per_shift[i=1:n], per_f(2 * pi * tspan[i] + u[8] * pi))
    @expression(model, per[i=1:n], per_f(2 * pi * tspan[i]))

    #Inlet Expressions
    #xin = @NLexpression(model, [c = 1:10, i = 1:n], 0.0)

    ts = LinRange(0, 2 * pi, 10000001)
    ratio = sum(rplug.methods.signal.(ts)) / length(ts)
    @expression(model, ndotin, rplug.properties.ndotin / (1 + ratio * u[4]))

    @expression(model, xin1[i=1:n], 0)
    @expression(model, xin2[i=1:n], u[1] * (1 - u[13] * per[i]))
    @expression(model, xin3[i=1:n], u[2] * (1 + u[5] * per[i]))
    @expression(model, xin4[i=1:n], u[3] * (1 - u[14] * per[i]))
    @expression(model, xin5[i=1:n], 0)
    @expression(model, xin6[i=1:n], u[9] * (1 - u[6] * per[i]))
    #@expression(model, xin7[i=1:n], x[7, i]) #Egal für optimierung
    @expression(model, xin8[i=1:n], u[10] * (1 + u[15] * per[i]))
    @expression(model, xin9[i=1:n], ndotin * (1 + u[4] * per_shift[i]))
    @expression(model, xin10[i=1:n], u[11] * (1 + u[15] * per[i]))

    #reflux
    if rplug.properties.reflux != 0
        nup = @expression(model, [i = 1:n], x[2, i, end] * x[9, i, end] + x[3, i, end] * x[9, i, end] + x[4, i, end] * x[9, i, end] + x[6, i, end] * x[9, i, end])
        xrec2 = @expression(model, [i = 1:n], x[2, i, end] * x[9, i, end] / nup[i])
        xrec3 = @expression(model, [i = 1:n], x[3, i, end] * x[9, i, end] / nup[i])
        xrec4 = @expression(model, [i = 1:n], x[4, i, end] * x[9, i, end] / nup[i])
        xrec6 = @expression(model, [i = 1:n], x[6, i, end] * x[9, i, end] / nup[i])
        nrec = @expression(model, [i = 1:n], rplug.properties.reflux * nup[i])
    end
    ## PFTR-Equations
    for i = 1:n
        for j = 1:m
            if i == 1
                for c = 1:9
                    @constraint(model, x[c, i, j] == x[c, end, j])
                end
            elseif j == 1
                if rplug.properties.reflux != 0
                    @NLconstraint(model, x[1, i, j] == xin1[i])
                    @NLconstraint(model, x[2, i, j] == (xin2[i] * xin9[i] + xrec2[i] * nrec[i]) / (xin9[i] + nrec[i]))
                    @NLconstraint(model, x[3, i, j] == (xin3[i] * xin9[i] + xrec3[i] * nrec[i]) / (xin9[i] + nrec[i]))
                    @NLconstraint(model, x[4, i, j] == (xin4[i] * xin9[i] + xrec4[i] * nrec[i]) / (xin9[i] + nrec[i]))
                    @NLconstraint(model, x[5, i, j] == xin5[i])
                    @NLconstraint(model, x[6, i, j] == (xin6[i] * xin9[i] + xrec6[i] * nrec[i]) / (xin9[i] + nrec[i]))
                    @NLconstraint(model, x[7, i, j] == 0)
                    @NLconstraint(model, x[8, i, j] == xin8[i])
                    @NLconstraint(model, x[9, i, j] == xin9[i] + nrec[i])
                else
                    xin = [xin1[i], xin2[i], xin3[i], xin4[i], xin5[i], xin6[i], 0, xin8[i], xin9[i], xin10[i]]
                    for c = 1:9
                        @NLconstraint(model, x[c, i, j] == xin[c])
                    end
                end
            else
                setRplugCons(model, dxt[:, i, j], dxz[:, i, j], x[:, i, j], xin10[i], rplug::Reactor)
            end
        end
    end

    ## Imlet Constraints
    @expression(model, mean_N2_in, sum((xin6[j] * xin9[j] * 0.5 + xin6[j-1] * xin9[j-1] * 0.5) * (tspan[j] - tspan[j-1]) for j = 2:n))
    @expression(model, mean_H2_in, sum((xin4[j] * xin9[j] * 0.5 + xin4[j-1] * xin9[j-1] * 0.5) * (tspan[j] - tspan[j-1]) for j = 2:n))

    setInputCons(model, u, mean_N2_in, u0, rplug, rplug.methods.optimize)
    n_in = xin9

    # for i = 1:n
    #     @NLconstraint(model, n_in[i] >= rplug.properties.ndotin * 0.1)
    # end
    ## Objective Functions
    n_in = xin9
    y_in2 = xin2
    y_in3 = xin3
    n_in2 = @expression(model, sum((n_in[j] * y_in2[j] + n_in[j-1] * y_in2[j-1]) * 0.5 * (tspan[j] - tspan[j-1]) for j = 2:n))
    n_in3 = @expression(model, sum((n_in[j] * y_in3[j] + n_in[j-1] * y_in3[j-1]) * 0.5 * (tspan[j] - tspan[j-1]) for j = 2:n))

    rplug.methods.obj1(model, x[:, :, end], tspan, n, rplug.properties.mcat)
    rplug.methods.obj2(model, x[:, :, end], n_in2, n_in3, tspan, n, rplug.properties.mcat)
    obj1 = model[:obj1]
    obj2 = model[:obj2]

    if rplug.methods.setObjective == 1
        @NLobjective(model, Max, obj1)
        @NLconstraint(model, obj2 >= epsilon)
    elseif rplug.methods.setObjective == 2
        @NLobjective(model, Max, obj2)
        @NLconstraint(model, obj2 >= epsilon)
    end
    # Solve
    @time solu = JuMP.optimize!(model)
    #rplug.methods.sol(model, iter)
    if rplug.methods.setObjective == 1
        if (termination_status(model) == TerminationStatusCode(4) || termination_status(model) == TerminationStatusCode(10) || termination_status(model) == TerminationStatusCode(19)) && value(model[:obj1]) > rplug.methods.sol.obj1[iter]
            rplug.methods.sol(model, iter)
        end
    else
        if (termination_status(model) == TerminationStatusCode(4) || termination_status(model) == TerminationStatusCode(10) || termination_status(model) == TerminationStatusCode(19)) && value(model[:obj2]) > rplug.methods.sol.obj2[iter]
            rplug.methods.sol(model, iter)
        end
    end
    return model
    MOI.Utilities.reset_optimizer(model)
    model = nothing
    GC.gc()
    nothing
end

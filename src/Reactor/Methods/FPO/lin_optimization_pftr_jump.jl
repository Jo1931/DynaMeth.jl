function opti_pftr_1(epsilon, param, setup)
    @unpack setop, ss, fixnin, start, sca, mold = setup
    @unpack n0, ordn, m0, ordm, nin, N2, per_f, u_lb, u_ub, l = param

    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 60.0 * 60, "max_iter" => 1000,
            "nlp_scaling_obj_target_gradient" => sca, "linear_solver" => "ma97",
            #"tol" => 1e-12,
            #"bound_relax_factor" => 0.0
        ), add_bridges=false)#

    n = n0 + (ordn - 2) * (n0 - 1)
    m = m0 + (ordm - 2) * (m0 - 1)
    ti = LinRange(0, 1, n0)
    zi = LinRange(0, 1, m0)
    nc = 3
    if start == true
        @variable(model, 0 <= x[c=1:nc, i=1:n, j=1:m] <= [1, 1, 1, 1, nin * 100][c], start = value.(mold[:x][c, i, j]))
        @variable(model, u_lb[i] <= u[i=1:8] <= u_ub[i], start = value.(mold[:u][i]))
    else
        @variable(model, 0 <= x[c=1:nc, i=1:n, j=1:m] <= [1, 1, 1, 1, nin * 100][c], start = 1)
        @variable(model, u_lb[i] <= u[i=1:8] <= u_ub[i])
    end
    tspa, zspa = GL_nodes(ti, zi, ordn, ordm)
    @expression(model, tspan[i=1:n], tspa[i])
    @expression(model, tend, u[6] * 3600)
    @expression(model, zspan[i=1:m], zspa[i])
    @expression(model, zend, l)

    dxt = dxdt(model, x, tspan, tend, ordn, n, m, nc)
    dxz = dxdz(model, x, zspan, zend, ordm, n, m, nc)

    register(model, :per_f, 1, per_f; autodiff=true)
    @expression(model, per_shift[i=1:n], per_f(2 * pi * tspan[i] + u[7] * 1 * pi))
    @expression(model, per[i=1:n], per_f(2 * pi * tspan[i]))


    @expression(model, xin1[i=1:n], u[1] * (1 + u[3] * per[i]))
    @expression(model, xin2[i=1:n], u[2] * (1 - u[4] * per[i]))
    @expression(model, xin3[i=1:n], 0)
    @expression(model, ndotin[i=1:n], u[8] * (1 + u[5] * per_shift[i]))



    @expression(model, mean_N2_in, sum((ndotin[j] * xin2[j] * 0.5 + ndotin[j-1] * xin2[j-1] * 0.5) * (tspan[j] - tspan[j-1]) for j = 2:n))
    @NLconstraint(model, abs(u[8] * N2 - mean_N2_in) <= 1e-5)

    @constraint(model, u[3] * u[1] - u[4] * u[2] == 0)
    @constraint(model, u[1] + u[2] == 1)
    if ss == true
        @constraint(model, u[3] == 0)
        @constraint(model, u[5] == 0)
    end
    if fixnin == true
        @constraint(model, u[8] == nin)
    end

    # @N



    ## PFTR-Equations
    for i = 1:n
        for j = 1:m
            if i == 1
                for c = 1:nc
                    @constraint(model, x[c, 1, j] == x[c, end, j])
                end
            elseif j == 1
                @constraint(model, x[1, i, 1] == xin1[i])
                @constraint(model, x[2, i, 1] == xin2[i])
                @constraint(model, x[3, i, 1] == xin3[i])
            else
                pftrcons1(model, dxt[:, i, j], dxz[:, i, j], x[:, i, j], ndotin[i], l)
            end
        end
    end
    n_in1 = @expression(model, sum((ndotin[j] * xin1[j] + ndotin[j-1] * xin1[j-1]) * 0.5 * (tspan[j] - tspan[j-1]) for j = 2:n))

    @expression(model, obj1, sum((ndotin[j] * x[3, j, end] * 0.5 + ndotin[j] * x[3, j-1, end] * 0.5) * (tspan[j] - tspan[j-1]) for j = 2:n))
    @expression(model, obj2, sum((ndotin[j] * x[3, j, end] * 0.5 + ndotin[j] * x[3, j-1, end] * 0.5) * (tspan[j] - tspan[j-1]) for j = 2:n) / n_in1)
    #@expression(model, obj2, sum(x[3, j]/(xin1[j]+1e-6) for j = 1:n)/n)
    if setop == 1
        @NLobjective(model, Max, obj1)
    else
        @NLobjective(model, Max, obj2)
    end
    @NLconstraint(model, obj2 >= epsilon)
    optimize!(model)
    return model
end
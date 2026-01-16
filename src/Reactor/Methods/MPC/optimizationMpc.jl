
sigmoid(x) = (atan(10000000 * x) + pi / 2) / pi#1/(1+exp(50*x))

function solve_opt_crl(cstr, y0, x0, u0, s0, uprev, tstart; sca=1.0, contr=0.0)

    te = cstr.methods.horizon
    ti = LinRange(0, 1, cstr.methods.n)
    tspan = GL_nodes(ti, cstr.methods.ordn)
    x_lb = [0, 0, 0, 0.05, 0, 0, 0, 0]
    x_ub = [1, 1, 1, 1, 1, 1, 0.95, cstr.properties.ndotin * 5]
    u_lb = [6.3177e-5, 0.00014876]
    u_ub = [0.00126234, 0.00198342]

    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 60.0 * 60, "max_iter" => 500,
            "nlp_scaling_obj_target_gradient" => sca,
            #"tol" => 1e-12,
            #"bound_relax_factor" => 0.0
        ), add_bridges=false)#,"acceptable_constr_viol_tol" => 1e-18))
    ## define Variables
    set_string_names_on_creation(model, false)

    @variable(model, x_lb[c] <= x[c=1:8, i=1:cstr.methods.nall] <= x_ub[c], start = x0[c, i])
    @variable(model, u_lb[c] <= u[c=1:2, i=1:cstr.methods.nall] <= u_ub[c], start = u0[c, i])
    @variable(model, u_lb[c] <= uset[c=1:2, i=1:cstr.methods.nall] <= u_ub[c], start = u0[c, i])
    @variable(model, z[i=1:cstr.methods.nall], start = 0)
    @variable(model, -1 <= s1[i=1:cstr.methods.nall] <= 1, start = s0[i])

    ## define grid (Gauß Lobatto points)


    # define period and frequency
    #@expression(model, tend, te)


    ## define derivative


    #return res
    dxt = dxdt(model, x, tspan, te, cstr.methods.ordn, cstr.methods.nall)

    ## periodic input function
    register(model, :sigmoid, 1, sigmoid; autodiff=true)


    nin = cstr.properties.ndotin
    yin = cstr.properties.yin

    #Inlet Expressions
    ndt = 0
    H2_dist = cstr.methods.disturb(tspan * te .+ tstart .- cstr.methods.deadtime)
    @expression(model, ts[i=1:cstr.methods.nall], tstart + tspan[i] * te)
    @expression(model, nin1[i=1:cstr.methods.nall], nin * yin[1])
    nin2 = Array{NonlinearExpr}(undef, cstr.methods.nall)
    nin3 = Array{NonlinearExpr}(undef, cstr.methods.nall)

    nin2[1:cstr.methods.nall] = @expression(model, [i = 1:cstr.methods.nall], u[1, i])
    nin3[1:cstr.methods.nall] = @expression(model, [i = 1:cstr.methods.nall], u[2, i])

    #@expression(model, nin3[i=1:cstr.methods.nall], u[2, i])
    #@expression(model, nin2[i=1:n], nin * yin[2])
    #@expression(model, nin3[i=1:n], nin * yin[3])
    @expression(model, nin4[i=1:cstr.methods.nall], nin * yin[4] * (1 .+ H2_dist[i]))
    @expression(model, nin5[i=1:cstr.methods.nall], nin * yin[5])
    @expression(model, nin6[i=1:cstr.methods.nall], nin * yin[6])

    ## CSTR-Equations


    [cstr.properties.ndotin * 0.3, cstr.properties.ndotin * 0.3]
    @constraint(model, u[1, 1] == cstr.methods.input.co2(tspan[1] * te .+ tstart))#uprev[1][end])
    @constraint(model, u[2, 1] == cstr.methods.input.co(tspan[1] * te .+ tstart))#uprev[2][end])  #@constraint(model, uset[1, 1] == uprev[1][end] * (1 + contr))
    #@constraint(model, uset[2, 1] == uprev[2][end] * (1 + contr))
    #@NLconstraint(model, 20 * (u[1, 1] - uprev[1]) / (tspan[2] - tspan[1]) / te == -uprev[1] + uset[1, 1])
    #@NLconstraint(model, 20 * (u[2, 1] - uprev[2]) / (tspan[2] - tspan[1]) / te == -uprev[2] + uset[1, 1])
    for c = 1:8
        @constraint(model, x[c, 1] == y0[c])
        #@constraint(model, 1.0 == yin[1] + yin[4] + yin[5] + yin[6] + u[1, 1] + u[2, 1])
        @constraint(model, z[1] == cstr.methods.err_int)
    end
    #@constraint(model, (u[1, 2] - u[1, 1]) * 1e5 == (uprev[1][end] - uprev[1][end-1]))
    #@constraint(model, (u[2, 2] - u[2, 1]) * 1e5 == (uprev[2][end] - uprev[2][end-1]))
    #@constraint(model, uset[1, 1] == uprev[1][end] * (1 + contr))
    #@constraint(model, uset[2, 1] == uprev[2][end] * (1 + contr))
    d = round(Int, cstr.methods.deadtime / (tspan[2] - tspan[1]) / te)  # Verzögerung in Schritten

    nin2_delay = Array{NonlinearExpr}(undef, cstr.methods.nall)
    nin3_delay = Array{NonlinearExpr}(undef, cstr.methods.nall)

    nin2_delay[1:d] = @expression(model, [i = 1:d], cstr.methods.input.co2(tspan[i] * te .+ tstart .- cstr.methods.deadtime))
    nin3_delay[1:d] = @expression(model, [i = 1:d], cstr.methods.input.co(tspan[i] * te .+ tstart .- cstr.methods.deadtime))

    nin2_delay[d+1:cstr.methods.nall] = @expression(model, [i = 1:cstr.methods.nall-d], u[1, i])
    nin3_delay[d+1:cstr.methods.nall] = @expression(model, [i = 1:cstr.methods.nall-d], u[2, i])
    T = 50
    @expression(model, nin8[i=1:cstr.methods.nall], nin1[i] + nin2_delay[i] + nin3_delay[i] + nin4[i] + nin5[i] + nin6[i])

    @NLexpression(model, obj2[i=1:cstr.methods.nall], x[8, i] * x[1, i] / (nin2_delay[i] + nin3_delay[i]))
    #@NLexpression(model, obj2[i=1:cstr.methods.nall], x[8, i] * x[1, i] / (nin2[i] + nin3[i]))

    #@NLexpression(model, obj2[i=1:cstr.methods.nall], x[1, i])

    @NLexpression(model, obj1[j=1:cstr.methods.nall], (x[8, j] * x[1, j]) * 6e4)

    for i = 2:cstr.methods.nall
        @NLconstraint(model, T * (u[1, i] - u[1, i-1]) / (tspan[2] - tspan[1]) / te == -u[1, i] + uset[1, i])
        @NLconstraint(model, T * (u[2, i] - u[2, i-1]) / (tspan[2] - tspan[1]) / te == -u[2, i] + uset[2, i])
    end
    for i = 2:cstr.methods.nall
        nin0 = [nin1[i], nin2_delay[i], nin3_delay[i], nin4[i], nin5[i], nin6[i], 0, nin8[i]]
        setCstrConsMpc(model, dxt[:, i], x[:, i], nin0, cstr::Reactor)
    end
    for i = 2:cstr.methods.nall
        @NLconstraint(model, obj2[i] + s1[i] == cstr.methods.setpoint)
        #@NLconstraint(model, obj2[i] + s1[i] == cstr.methods.setpoint)
        #@NLconstraint(model, obj1[i] + s1[i] == 0.9)
    end

    for i = 2:cstr.methods.nall
        @NLconstraint(model, 500 * (z[i] - z[i-1]) / (tspan[2] - tspan[1]) / te == (obj2[i] - cstr.methods.setpoint))
    end

    #@expression(model, penalty, -sum((u[1, i] / nin - u[1, i-1] / nin)^2 + (u[2, i] / nin - u[2, i-1] / nin)^2 for i in 2:cstr.methods.nall))

    #@expression(model, obj1, sum((x[8, j] * x[1, j] + x[8, j-1] * x[1, j-1]) * (tspan[j] - tspan[j-1]) for j = 2:cstr.methods.nall) * 6e4)

    #@NLobjective(model, Max, 1 * sum(obj1[i] for i = 1:cstr.methods.nall) + 0 * penalty - 1000000 * sum(sigmoid(s1[i]) * s1[i]^2 for i in 1:cstr.methods.nall) - 0 * sum(z[i]^2 for i in 1:cstr.methods.nall))
    #@NLobjective(model, Max, -sum(4 * z[i]^2 + 0 * (obj2[i] - cstr.methods.setpoint)^2 for i in 1:cstr.methods.nall))
    if set_object == 1
        #@NLobjective(model, Max, 1 * sum(obj1[i] for i = 1:cstr.methods.nall) - 1000000 * sum(sigmoid(s1[i]) * s1[i]^2 for i in 1:cstr.methods.nall))
        @NLobjective(model, Max, 1 * sum(obj1[i] for i = d+1:cstr.methods.nall) - 1e6 * sum(s1[i]^2 for i in d+1:cstr.methods.nall))

    elseif set_object == 2
        @NLobjective(model, Max, -1e5 * sum((obj2[i] - cstr.methods.setpoint)^2 for i in 3:cstr.methods.nall))
    elseif set_object == 3
        @NLobjective(model, Max, -1e3 * sum(z[i] for i in 1:cstr.methods.nall)^2 - 1 * 1e5 * sum(s1[i]^2 for i in 1:cstr.methods.nall) + 0 * sum(obj1[i] for i = 1:cstr.methods.nall))
        #@NLobjective(model, Max, 1 * sum(obj1[i] for i = 1:cstr.methods.nall) - 1000000 * sum(4 * z[i]^2 + 1 * sigmoid(s1[i]) * s1[i]^2 for i in 1:cstr.methods.nall))

    end
    @time solu = JuMP.optimize!(model)

    return model
end
function simulate_plant(cstr_plant, y0, ts, saveat=LinRange(ts[1], ts[end], 2))
    dy0 = zero(y0)   # Startwerte für die Ableitungen (häufig Null)
    differential_vars = [true, true, true, true, true, true, true, false]
    prob = DAEProblem(dae_system_mpc, dy0, y0, ts, cstr_plant, differential_vars=differential_vars, saveat=saveat)
    sol = solve(prob, DFBDF())
    x = sol[:, :] .* (0.00 * (2 * rand(length(sol[:, 1]), length(sol[1, :])) .- 1) .+ 1)
    Y = x[1, :] .* x[8, :] ./ (cstr_plant.methods.input.co(sol.t .- cstr_plant.methods.deadtime) + cstr_plant.methods.input.co2(sol.t .- cstr_plant.methods.deadtime))
    return x, Y, sol.t, sol[:, end]
end

function objective_function(x, p)
    fac = p.fac
    solMPC = p.solMPC
    cstr_est = p.cstr
    cstr_est.methods.para = ones(3) * x
    lp = 2
    len = length(solMPC.t)
    if len <= lp
        n = len - 1
    else
        n = lp - 1
    end
    tsave = solMPC.t[end-n:end]

    y0 = solMPC.x[:, end-n]
    dy0 = zero(y0)   # Startwerte für die Ableitungen (häufig Null)
    differential_vars = [true, true, true, true, true, true, true, false]
    saveat = tsave
    prob = DAEProblem(dae_system_mpc, dy0, y0, (tsave[1], tsave[end]), cstr_est, differential_vars=differential_vars, saveat=saveat)
    sol = solve(prob, DFBDF())
    xx = sol[:, :]
    #Y = xx[1, :] .* xx[8, :] ./ (cstr_est.methods.input.co(sol.t) + cstr_est.methods.input.co2(sol.t))
    YMess = solMPC.obj2[end-n:end]
    XMess = solMPC.x[1:3, end-n:end]
    #x = sol[:, :]
    #Y = x[1, :] .* x[8, :] ./ (cstr_plant.methods.input.co(sol.t .- cstr_plant.methods.deadtime) + cstr_plant.methods.input.co2(sol.t .- cstr_plant.methods.deadtime))

    #sum((Y .- YMess) .^ 2)
    #print(" ")
    #print(sum(((X .- XMess) ./ (XMess .+ 0.1)) .^ 2))
    #print(" ")
    sum(((sol[1:3, :] .- XMess)) .^ 2)

end



function estimate_p(x0, fac, solMPC, cstr_est)
    p = (fac=fac, solMPC=solMPC, cstr=cstr_est)
    lower = 0.2
    upper = 5.0

    f(x) = objective_function(x, p)
    sol = optimize(f, lower, upper, Optim.GoldenSection())
    return ones(3) * Optim.minimizer(sol)
end

function objective_function_2(x, p)
    fac = p.fac
    solMPC = p.solMPC
    cstr_est = p.cstr
    cstr_est.methods.para = ones(3)
    cstr_est.methods.para[1] = x[1]
    cstr_est.methods.para[2] = x[2]
    cstr_est.methods.para[3] = x[3]
    lp = 3
    len = length(solMPC.t)
    if len <= lp
        n = len - 1
    else
        n = lp - 1
    end
    tsave = solMPC.t[end-n:end]

    y0 = solMPC.x[:, end-n]
    dy0 = zero(y0)   # Startwerte für die Ableitungen (häufig Null)
    differential_vars = [true, true, true, true, true, true, true, false]
    saveat = tsave
    prob = DAEProblem(dae_system_mpc, dy0, y0, (tsave[1], tsave[end]), cstr_est, differential_vars=differential_vars, saveat=saveat)
    sol = solve(prob, DFBDF())
    xx = sol[:, :]
    #Y = xx[1, :] .* xx[8, :] ./ (cstr_est.methods.input.co(sol.t) + cstr_est.methods.input.co2(sol.t))
    YMess = solMPC.obj2[end-n:end]
    XMess = solMPC.x[[1, 2, 5], end-n:end]
    xx = sol[:, :]
    Y = xx[1, :] .* xx[8, :] ./ (cstr_plant.methods.input.co(sol.t .- cstr_plant.methods.deadtime) + cstr_plant.methods.input.co2(sol.t .- cstr_plant.methods.deadtime))

    #1e5 * sum((Y .- YMess) .^ 2) + fac * ((1 - x[1])^2 + (1 - x[3])^2)
    #print(" ")
    #print(sum(((X .- XMess) ./ (XMess .+ 0.1)) .^ 2))
    #print(" ")
    1e7 * sum(((sol[[1, 2, 5], :] .- XMess)) .^ 2) + fac * ((1 - x[1])^2 + (1 - x[3])^2)

end

using Optim

function estimate_p_2(x0, fac, solMPC, cstr_est)
    p = (fac=fac, solMPC=solMPC, cstr=cstr_est)

    # untere und obere Schranken für 2 Parameter
    lower = [0.01, 0.01, 0.01]
    upper = [20.0, 20.0, 20.0]

    # Zielfunktion
    f(x) = objective_function_2(x, p)
    opt = Optim.Options(time_limit=10.0, show_trace=true)
    # Optimierung mit Fminbox + Nelder–Mead
    sol = optimize(f, lower, upper, x0,
        Fminbox(LBFGS()), opt)

    return Optim.minimizer(sol)[1:3]
end

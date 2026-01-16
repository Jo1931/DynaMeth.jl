function solve_opt(cstr; solver=Optim.NelderMead(), x0=ComponentArray(get_kinetic_param(cstr.kinetic.parameter))[1:end] .+ 1e-20, iter=1000, show_trace=true)
    #x0_unwrapped = ComponentArray(get_kinetic_param(cstr.kinetic.parameter))[1:end]
    x0_ = wrappedParam(x0, cstr)

    optf = OptimizationFunction(cstr.methods.objectiveSS)
    prob = Optimization.OptimizationProblem(optf, x0_, cstr, show_trace=show_trace, iterations=iter)
    @time sol = solve(prob, solver)
    #@time sol = solve(prob, Optim.ParticleSwarm(; n_particles=20), show_trace=true, iterations=10000)

    solu = unwrappedParam(sol.u, cstr)
    set_kinetic_param!(cstr, solu)
    return solu
end
function obnomad(x, cstr)
    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    u = unwrappedParam(x, cstr)
    set_kinetic_param!(cstr, u)
    xsim = solveSS(cstr, i)
    xmess = cstr.methods.xout_mess[i]
    out = sum(((xsim[1:3] - xmess[1:3]) ./ (xmess[1:3] .+ 1e-1)) .^ 2)
    #out = sum(((xsim[1:3] - xmess[1:3]) .^ 2))
    for i in experimental_set[2:end]
        xsims = solveSS(cstr, i)
        xmesss = cstr.methods.xout_mess[i]
        out = out + sum(((xsims[1:3] - xmesss[1:3]) ./ (xmesss[1:3] .+ 1e-1)) .^ 2)
        #out = out + sum(((xsims[1:3] - xmesss[1:3]) .^ 2))
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end

    #cstr.methods.xout_sim = xsim
    res = out
end
function obnomad_dyn(x, cstr)
    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    u = unwrappedParam(x, cstr)
    set_kinetic_param!(cstr, u)
    xsim = solveSS(cstr, i)
    xmess = cstr.methods.xout_mess[i]
    out = sum(((xsim[1:3] - xmess[1:3]) ./ (xmess[1:3] .+ 1e-1)) .^ 2)
    #out = sum(((xsim[1:3] - xmess[1:3]) .^ 2))
    for i in experimental_set[2:end]
        xsims = solveSS(cstr, i)
        xmesss = cstr.methods.xout_mess[i]
        out = out + sum(((xsims[1:3] - xmesss[1:3]) ./ (xmesss[1:3] .+ 1e-1)) .^ 2)
        #out = out + sum(((xsims[1:3] - xmesss[1:3]) .^ 2))
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end
    data1 = cstr.methods.files

    ts = data1[!, "dynamic_1"] * 60
    ts, x_out, sol, xin = simulateVollbrechtDyn(cstr, ts)

    ch3oh = ((sol[1, :] - data1[!, "dynamic_2"]) ./ (data1[!, "dynamic_2"] .+ 1e-1)) .^ 2
    co2 = ((sol[2, :] - data1[!, "dynamic_3"]) ./ (data1[!, "dynamic_3"] .+ 1e-1)) .^ 2
    co = ((sol[3, :] - data1[!, "dynamic_4"]) ./ (data1[!, "dynamic_3"] .+ 1e-1)) .^ 2
    dyn = sum(ch3oh .^ 2) + sum(co2 .^ 3) + sum(co .^ 4)

    out = out .+ dyn * 100
    #cstr.methods.xout_sim = xsim
    res = out + regularization(u)
end
function solve_opt_NOMAD(cstr; x0=ComponentArray(get_kinetic_param(cstr.kinetic.parameter))[1:end] .+ 1e-20, show_trace=true)
    #x0_unwrapped = ComponentArray(get_kinetic_param(cstr.kinetic.parameter))[1:end]
    x0_ = wrappedParam(x0, cstr)

    optf = OptimizationFunction(obnomad)
    prob = Optimization.OptimizationProblem(optf, x0_, cstr, show_trace=true, iterations=10000)
    sol = solve(prob, NOMADOpt())
    solu = unwrappedParam(sol.u, cstr)
    set_kinetic_param!(cstr, solu)
    return solu
end
function solve_opt_Rexp(cstr; solver=Optim.NelderMead(), x0=ComponentArray(get_kinetic_param(cstr.kinetic.parameter))[1:end], iter=1000, show_trace=true, diff=AutoForwardDiff())
    #x0_unwrapped = ComponentArray(get_kinetic_param(cstr.kinetic.parameter))[1:end]
    x0_ = wrappedParam(x0, cstr)

    optf = OptimizationFunction(obj_Rexp, diff)
    prob = Optimization.OptimizationProblem(optf, x0_, cstr, show_trace=show_trace, iterations=iter)
    @time sol = solve(prob, solver)
    #@time sol = solve(prob, Optim.ParticleSwarm(; n_particles=20), show_trace=true, iterations=10000)

    solu = unwrappedParam(sol.u, cstr)
    set_kinetic_param!(cstr, solu)
    return solu
end


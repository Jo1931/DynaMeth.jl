function simulateExperiments(p::Reactor, file)
    xin, x_out, t, t0 = readFile(file)
    idx = p.experiments.idx
    spl_CO2 = Spline1D(t, xin[2], k=1)
    spl_N2 = Spline1D(t, xin[6], k=1)
    spl_H2 = Spline1D(t, xin[4], k=1)
    spl_CO = Spline1D(t, xin[3], k=1)
    spl_n = Spline1D(t, xin[7], k=1)
    spl_T = Spline1D(t, xin[8], k=1)
    spl_P = Spline1D(t, xin[9], k=1)

    spl = [spl_CO2, spl_CO, spl_H2, spl_N2, spl_n, spl_T, spl_P]
    x0 = [x_out["CH3OH"][idx], x_out["CO2"][idx], x_out["CO"][idx], x_out["H2"][idx], x_out["H2O"][idx], x_out["N2"][idx], 0.7]


    Solver = DFBDF()
    #Solver=IDA()
    #Solver= Rodas5P()
    tend = t[end]
    u0 = x0
    du0 = zeros(7)
    differential_vars = [true, true, true, true, true, true, true] # like a mass-matrix
    tspan = (t0[idx], tend)
    ts = t0[idx:end, 1]
    p.properties.activity = xin[10]
    pn = (spl, p)
    prob = DAEProblem(daeFuncMess, du0, u0, tspan, differential_vars=differential_vars) #defin Problem
    sol = try
        DifferentialEquations.solve(prob, Solver, p=pn, saveat=ts, maxiters=5e3) # Solve Problem        
    catch
        x = zeros(7, length(ts))
        for i = 1:7
            x[i, :] = ones(1, length(ts)) * u0[i]
        end
        x
    end
    if length(sol[1, :]) < length(ts)
        x = zeros(7, length(ts))
        for i = 1:7
            x[i, :] = ones(1, length(ts)) * u0[i]
        end
        x
        sol = x
    end
    return ts, x_out, sol, xin
end
function simulateVollbrechtDyn(p::Reactor, ts=LinRange(0, 350 * 60, 500))
    data = p.methods.files
    T = 250 + 273.15
    P = 50 * 1e5
    n = (240 / 60 * 1e-6 * p.properties.pN) / (8.314 * p.properties.TN)
    t = [0.0, 140 * 60, 140 * 60 + 1, 210 * 60, 210 * 60 + 1, 280 * 60, 280 * 60 + 1, 350 * 60]
    xin2 = [0.0, 0.0, 0.119, 0.119, 0.0, 0.0, 0.119, 0.119]
    xin3 = [0.126, 0.126, 0.0, 0.0, 0.126, 0.126, 0.0, 0.0]
    xin4 = [0.722, 0.722, 0.715, 0.715, 0.722, 0.722, 0.715, 0.715]
    xin6 = [0.16, 0.16, 0.166, 0.166, 0.16, 0.16, 0.166, 0.166]
    xin7 = [n, n, n, n, n, n, n, n]
    xin8 = [T, T, T, T, T, T, T, T]
    xin9 = [P, P, P, P, P, P, P, P]
    spl_CO2 = Spline1D(t, xin2, k=1)
    spl_N2 = Spline1D(t, xin6, k=1)
    spl_H2 = Spline1D(t, xin4, k=1)
    spl_CO = Spline1D(t, xin3, k=1)
    spl_n = Spline1D(t, xin7, k=1)
    spl_T = Spline1D(t, xin8, k=1)
    spl_P = Spline1D(t, xin9, k=1)

    spl = [spl_CO2, spl_CO, spl_H2, spl_N2, spl_n, spl_T, spl_P]
    ##todo Vollbrecht daten Einlese
    x0 = [data[!, "dynamic_2"][1], data[!, "dynamic_3"][1], data[!, "dynamic_4"][1], data[!, "dynamic_5"][1], data[!, "dynamic_6"][1], data[!, "dynamic_7"][1], 0.9]


    Solver = DFBDF()
    #Solver = IDA()
    tend = t[end]
    u0 = x0
    du0 = zeros(7)
    differential_vars = [true, true, true, true, true, true, true] # like a mass-matrix
    tspan = (data[!, "dynamic_1"][1], data[!, "dynamic_1"][end] * 60)
    #ts = LinRange(tspan[1],tspan[2],200)
    pn = (spl, p)
    prob = DAEProblem(daeFuncMess, du0, u0, tspan, differential_vars=differential_vars) #defin Problem
    sol = try
        DifferentialEquations.solve(prob, Solver, p=pn, saveat=ts, maxiters=5e3, dtmax=100.0)
    catch
        x = zeros(7, length(ts))
        for i = 1:7
            x[i, :] = ones(1, length(ts)) * u0[i]
        end
        x * NaN
    end  # Solve Problem
    if length(sol[1, :]) < length(ts)
        x = zeros(7, length(ts))
        for i = 1:7
            x[i, :] = ones(1, length(ts)) * u0[i]
        end
        sol = x * NaN
    end
    return ts, data, sol, []
end

function simulateAllExperiments(p::Reactor, Files)
    idx = p.experiments.idx
    tv = []
    xoutv = []
    solv = []
    xinv = []
    err = []
    for file in Files
        t, x_out, sol, xin = simulateExperiments(p::Reactor, file)
        e1 = sum(sqrt.(((sol[1, :] - x_out["CH3OH"][idx:end]) ./ x_out["CH3OH"][idx:end]) .^ 2)) / length(t)
        e2 = sum(sqrt.(((sol[2, :] - x_out["CO2"][idx:end]) ./ x_out["CO2"][idx:end]) .^ 2)) / length(t)
        e3 = sum(sqrt.(((sol[3, :] - x_out["CO"][idx:end]) ./ x_out["CO"][idx:end]) .^ 2)) / length(t)
        e4 = sum(sqrt.(((sol[4, :] - x_out["H2"][idx:end]) ./ x_out["H2"][idx:end]) .^ 2)) / length(t)
        e5 = sum(sqrt.(((sol[5, :] - x_out["H2O"][idx:end]) ./ x_out["H2O"][idx:end]) .^ 2)) / length(t)
        e6 = sum(sqrt.(((sol[6, :] - x_out["N2"][idx:end]) ./ x_out["N2"][idx:end]) .^ 2)) / length(t)
        push!(tv, t)
        push!(xoutv, x_out)
        push!(solv, sol)
        push!(xinv, xin)
        push!(err, (e1 + e2 + e3) / 3)
        print(1)
    end
    return tv, xoutv, solv, xinv, err
end
function simulateExperimentsSteady(p::Reactor, file, xinn)
    xin, x_out, t, t0 = readFile(file)
    idx = p.experiments.idx
    spl_CO2 = Spline1D(t, xin[2], k=1)
    spl_N2 = Spline1D(t, xin[6], k=1)
    spl_H2 = Spline1D(t, xin[4], k=1)
    spl_CO = Spline1D(t, xin[3], k=1)
    spl_n = Spline1D(t, xin[7], k=1)
    spl_T = Spline1D(t, xin[8], k=1)
    spl_P = Spline1D(t, xin[9], k=1)

    spl = [spl_CO2, spl_CO, spl_H2, spl_N2, spl_n, spl_T, spl_P]
    x0 = [x_out["CH3OH"][idx], x_out["CO2"][idx], x_out["CO"][idx], x_out["H2"][idx], x_out["H2O"][idx], x_out["N2"][idx], 0.7]


    Solver = DFBDF()
    #Solver=IDA()
    tend = t[end]
    u0 = x0
    du0 = zeros(7)
    differential_vars = [true, true, true, true, true, true, true] # like a mass-matrix
    tspan = (t0[idx], tend)
    ts = t0[idx:end, 1]
    pn = (xinn, p)
    prob = DAEProblem(daeFuncMessSteady, du0, u0, tspan, differential_vars=differential_vars) #defin Problem

    sol = DifferentialEquations.solve(prob, Solver, p=pn, saveat=ts) # Solve Problem
    return ts, x_out, sol, xin
end
function daeFuncMessSteady(out, dy, y, reac, t)
    ## calc inlets
    RHS = calcRHSexp(y, reac)
    LHS = calcLHSexp(y, reac)
    out .= LHS * dy - RHS
end
function simulateSS(cstr::Reactor, idx)
    Solver = DFBDF()
    tend = 30000
    cstr.properties.yin = vcat(cstr.methods.xin_mess[idx][1:6], 0.1)
    cstr.properties.ndotin = cstr.methods.xin_mess[idx][7]
    cstr.properties.T = cstr.methods.xin_mess[idx][8]
    cstr.properties.P = cstr.methods.xin_mess[idx][9]
    cstr.properties.activity = cstr.methods.xin_mess[idx][10]


    u0 = cstr.properties.yin
    du0 = zeros(7)
    differential_vars = [true, true, true, true, true, true, true] # like a mass-matrix
    tspan = (0, tend)
    pn = cstr
    prob = DAEProblem(daeFuncMessSteady, du0, u0, tspan, differential_vars=differential_vars) #defin Problem
    sol = DifferentialEquations.solve(prob, Solver, p=pn)
    return sol[1:7, end]#, cstr.methods.xout_mess[idx]

end
function solveSS(cstr::Tr, idx; param=get_kinetic_param(cstr.kinetic.parameter)) where {Tr<:Union{Reactor,ComponentVector}}
    yin = vcat(cstr.methods.xin_mess[idx][1:6], 0.6)
    nin = cstr.methods.xin_mess[idx][7]
    T = cstr.methods.xin_mess[idx][8]
    P = cstr.methods.xin_mess[idx][9]
    activity = cstr.methods.xin_mess[idx][10]
    u0 = vcat(cstr.methods.xout_mess[idx][1:6], 0.6)

    f(y, p) = calcRHSexp(y, p, T, P * 1e-5, yin, nin, param, activity)
    prob = NonlinearProblem(f, u0, cstr)

    sol = solve(prob, NewtonRaphson())
    if sol.retcode != ReturnCode.Success
        sol = u0 * 0.01
    end
    #sol = solve(prob, NewtonRaphson())
    return sol#, cstr.methods.xout_mess[idx]

end
function solveSS(cstr::Reactor,)
    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    xsim = solveSS(cstr, i)
    xmess = cstr.methods.xout_mess[i]
    for i in experimental_set[2:end]
        xsims = solveSS(cstr, i)
        xmesss = cstr.methods.xout_mess[i]
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end
    return xsim, xmess
end
function simulateSS(cstr::Reactor,)
    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    xsim = simulateSS(cstr, i)
    xmess = cstr.methods.xout_mess[i]
    for i in experimental_set[2:end]
        xsims = simulateSS(cstr, i)
        xmesss = cstr.methods.xout_mess[i]
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end
    return xsim, xmess
end
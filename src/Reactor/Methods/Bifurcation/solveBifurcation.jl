
function solveBifSS(cstr::Reactor, x0=vcat(cstr.properties.yin, 0.5, cstr.properties.Tc, cstr.properties.ndotin))
    f = calcRHSsingle
    prob = NonlinearProblem(f, x0, cstr)
    solve(prob, NewtonRaphson())
end
function solveBifSSCascade(cstr::Reactor, x0=repeat(vcat(cstr.properties.yin, 0.5, cstr.properties.Tc, cstr.properties.ndotin), cstr.properties.num_casc))
    f = calcRHScascade
    prob = NonlinearProblem(f, x0, cstr)
    solve(prob, NewtonRaphson())
end
function solveBifArc(cstr::Reactor)
    @unpack Tcprev, xprev, Tcpred, xpred, xpreprev, Tcpreprev, ds = cstr.methods.predictor
    x0 = vcat(xpred, Tcpred)
    f = calcRHSsingleArc
    prob = NonlinearProblem(f, x0, cstr)
    sol = solve(prob, NewtonRaphson())
    return sol.u[1:9], sol.u[10]
end
function solveBifArcCascade(cstr::Reactor)
    @unpack Tcprev, xprev, Tcpred, xpred, xpreprev, Tcpreprev, ds = cstr.methods.predictor
    x0 = vcat(xpred, Tcpred)
    f = calcRHSCascadeArc
    prob = NonlinearProblem(f, x0, cstr)
    sol = solve(prob, NewtonRaphson())#, abstol=1e-18, reltol=1e-18)
    return sol.u[1:9*cstr.properties.num_casc], sol.u[end]
end
function solveIsothermal(Tc, cstr::Reactor, x0=vcat(cstr.properties.yin, 0.5, cstr.properties.Tc, cstr.properties.ndotin))
    cstr.properties.Tc = Tc
    f = calcRHSbifIsothermal
    prob = NonlinearProblem(f, x0, cstr)
    solve(prob, NewtonRaphson())
end
function solveIsothermalCascade(Tc, cstr::Reactor, x0=vcat(cstr.properties.yin, 0.5, cstr.properties.Tc, cstr.properties.ndotin)' .* ones(cstr.properties.num_casc, 9))
    cstr.properties.Tc = Tc
    f = calcRHSCascadeIsothermal
    x00 = reshape(x0, length(x0))
    prob = NonlinearProblem(f, x00, cstr)
    solve(prob, NewtonRaphson())
end

function solveHeat(Tvec, Tc0, cstr::Reactor)
    solv = solveIsothermal(Tvec[1], cstr)
    Heat_prodv = calcHeatProduction(solv.u, Tvec[1], cstr.properties.mcat, cstr)
    Heat_remov = calcHeatRemoval(Tvec[1], Tc0, Tc0, cstr.properties.ndotin, cstr.properties.A, cstr)
    x0 = solv.u
    for Tc in Tvec[2:end]
        sol = solveIsothermal(Tc, cstr::Reactor, x0)
        Heat_prod = calcHeatProduction(sol, Tc, cstr.properties.mcat, cstr)
        Heat_remo = calcHeatRemoval(Tc, Tc0, Tc0, cstr.properties.ndotin, cstr.properties.A, cstr)
        x0 = sol.u
        solv = hcat(solv, sol)
        Heat_prodv = vcat(Heat_prodv, Heat_prod)
        Heat_remov = vcat(Heat_remov, Heat_remo)
    end
    return Heat_prodv, Heat_remov, solv
end

function solveHeatCascade(Tvec, Tc0, cstr::Reactor)
    num_casc = cstr.properties.num_casc
    #cstr.properties.Tc = Tvec[1]

    solv = solveIsothermalCascade(Tvec[1], cstr)
    count = 0
    while ReturnCode.Success != solv.retcode && count < 100
        x0 = solv.u
        x0[1:6] = rand(6)
        solv = solveIsothermalCascade(Tvec[1], cstr)
        count = count + 1
    end
    if count == 1000
        error("to many iterations")
    end
    xm = reshape(solv.u, num_casc, 9)

    Heat_prodv = zeros(length(Tvec), num_casc)
    Heat_remov = zeros(length(Tvec), num_casc)
    j = 1
    for i = 1:num_casc
        Heat_prodv[j, i] = calcHeatProduction(xm[i, :], Tvec[1], cstr.properties.mcat / num_casc, cstr)
        if i == 1
            ndotin = cstr.properties.ndotin
            Tin = Tc0
        elseif i == 2
            Tin = 540
            ndotin = xm[i-1, 9]
        else
            Tin = Tc0
            ndotin = xm[i-1, 9]
        end
        Heat_remov[j, i] = calcHeatRemoval(Tvec[1], Tin, Tc0, ndotin, cstr.properties.A / num_casc, cstr)
    end

    for Tc in Tvec[2:end]
        j = j + 1
        sol = solveIsothermalCascade(Tc, cstr::Reactor, xm)
        xm = reshape(sol.u, num_casc, 9)
        solv = hcat(solv, sol)
        for i = 1:num_casc
            Heat_prodv[j, i] = calcHeatProduction(xm[i, :], Tc, cstr.properties.mcat / num_casc, cstr)
            if i == 1
                ndotin = cstr.properties.ndotin
                Tin = Tc0
            elseif i == 2
                Tin = 540
                ndotin = xm[i-1, 9]
            else
                Tin = Tc0
                ndotin = xm[i-1, 9]
            end
            Heat_remov[j, i] = calcHeatRemoval(Tc, Tin, Tc0, ndotin, cstr.properties.A / num_casc, cstr)
        end
    end
    return Heat_prodv, Heat_remov, solv
end

function solveContinuationSingle(T0, Tend, cstr::Reactor; ds=1)

    cstr.properties.Tc = T0


    nin = cstr.properties.ndotin
    yin = cstr.properties.yin
    sol = solveIsothermal(T0, cstr)
    sol = solveBifSS(cstr, sol.u)
    y = sol.u
    sols = Solution_2d()


    J = calcSysMatrix(cstr, y)
    sols.eig = eigvals(J)
    push!(sols.y, y)
    push!(sols.obj1, y[9] * y[1])
    push!(sols.obj2, y[9] * y[1] / (nin * yin[2] + nin * yin[3]))
    push!(sols.Tv, T0)
    push!(sols.nin, nin)


    cstr.methods.predictor = Predictor(y, T0, ds)

    y, Tc = solveBifArc(cstr)
    cstr.methods.predictor(y, Tc)
    while Tc < Tend
        y, Tc = solveBifArc(cstr)
        cstr.properties.Tc = Tc
        cstr.methods.predictor(y, Tc)
        A = calcSysMatrix(cstr, y)

        sols.eig = hcat(sols.eig, eigvals(A))
        push!(sols.y, y)
        push!(sols.obj1, y[9] * y[1])
        push!(sols.obj2, y[9] * y[1] / (nin * yin[2] + nin * yin[3]))
        push!(sols.Tv, Tc)
        push!(sols.nin, nin)
    end
    return sols
end
function solveContinuationCascade(T0, Tend, cstr::Reactor; ds=1, smax=100000)

    cstr.properties.Tc = T0

    nin = cstr.properties.ndotin
    yin = cstr.properties.yin

    x0 = repeat(vcat(cstr.properties.yin, 0.5, cstr.properties.Tc, cstr.properties.ndotin), cstr.properties.num_casc)

    sol = solveBifSSCascade(cstr::Reactor, x0)
    ym = reshape(sol.u, 9, cstr.properties.num_casc)'

    y = sol.u

    sols = Solution_2d()


    J = calcSysMatrixCascade(cstr, y)
    sols.eig = eigvals(J)
    push!(sols.y, y)
    push!(sols.obj1, ym[end, 9] * ym[end, 1])
    push!(sols.obj2, ym[end, 9] * ym[end, 1] / (nin * yin[2] + nin * yin[3]))
    push!(sols.Tv, T0)
    push!(sols.nin, nin)


    cstr.methods.predictor = Predictor(y, T0, ds)

    y, Tc = solveBifArcCascade(cstr)
    cstr.methods.predictor(y, Tc)
    counter = 0
    while Tc < Tend && counter < smax
        counter = counter + 1
        if mod(counter, 100) == 0
            print(counter)
        end
        y, Tc = solveBifArcCascade(cstr)
        cstr.properties.Tc = Tc
        J = calcSysMatrixCascade(cstr, y)
        ym = reshape(y, 9, cstr.properties.num_casc)'
        cstr.methods.predictor(y, Tc)

        sols.eig = hcat(sols.eig, eigvals(J))
        push!(sols.y, y)
        push!(sols.obj1, ym[end, 9] * ym[end, 1])
        push!(sols.obj2, ym[end, 9] * ym[end, 1] / (nin * yin[2] + nin * yin[3]))

        push!(sols.Tv, Tc)
        push!(sols.nin, nin)
    end
    return sols
end


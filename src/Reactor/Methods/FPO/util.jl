function squarewave(x; period=1.0, sharpness=100.0)
    tanh.(sharpness * sin.(2π * x / period))
end
function collect_results_JuMP(model::Model)
    d = copy(model.obj_dict)
    for (k, v) in d
        d[k] = value.(v)
    end
    return d
end

function insert_result!(full::Dict, all::Dict, i::Int)
    testkey = :x
    if size(full[testkey]) == size(all[testkey])

        for k in keys(full)
            try
                dims = size(full[k])
                new = zeros((dims..., i))
                nd = ndims(new)
                inds = ntuple(_ -> Colon(), nd - 1)  # erzeugt (Colon(), Colon(), ...)
                new[inds..., 1] = full[k]
                new[inds..., i] = all[k]
                full[k] = new
            catch
            end
        end
    else
        for k in keys(full)
            try
                full[k][ntuple(_ -> Colon(), ndims(full[k]) - 1)..., i] = all[k]
            catch
            end
        end
    end
end
function yield_fbr(u, cstr)
    u[1] * u[9] / (cstr.properties.yin[2] + cstr.properties.yin[3]) / cstr.properties.ndotin
end

function calc_equi_yield(sol, P=5e6)
    i = 1
    bifurcation = SetupBifurcationAnalysis(kinetic=Seidel())
    cstr = Reactor(bifurcation)
    cstr.properties.P = P
    cstr.properties.yin = sol.xin[1:6, 1, i]
    cstr.properties.ndotin = 1e-1
    cstr.properties.Tc = sol.xin[8, 1, i]
    cstr.properties.Tin = sol.xin[8, 1, i]

    soln = solveIsothermal(200, cstr)
    #cstr.properties.ndotin = 1e-3
    soln = solveIsothermal(sol.xin[8, 1, i], cstr, soln.u)
    kk = LinRange(2, 10, 20)
    for j in kk
        cstr.properties.ndotin = 1 / 10^j
        soln = solveIsothermal(sol.xin[8, 1, i], cstr, soln.u)
    end

    yield_fbr(soln.u, cstr)
    Y = zeros(40)


    for i = 1:40
        cstr.properties.yin = sol.xin[1:6, 1, i]
        cstr.properties.Tc = sol.xin[8, 1, i]
        cstr.properties.Tin = sol.xin[8, 1, i]




        #sol=solveIsothermal(200, cstr)#, x0=vcat(cstr.properties.yin, 0.5, cstr.properties.Tc, cstr.properties.ndotin))

        soln = solveIsothermal(sol.xin[8, 1, i], cstr, soln.u)

        Y[i] = yield_fbr(soln.u, cstr)
    end

    return Y
end

function calc_equi_yield_dyn(sol, P=5e6)
    Ym = zeros(40, length(sol.xin[1, :, 1]))
    for j = 1:length(sol.xin[1, :, 1])
        i = 1
        bifurcation = SetupBifurcationAnalysis(kinetic=Seidel())
        cstr = Reactor(bifurcation)
        cstr.properties.P = P
        cstr.properties.yin = sol.xin[1:6, j, i]
        cstr.properties.ndotin = 1e-3
        cstr.properties.Tc = sol.xin[8, j, i]


        soln = solveIsothermal(200, cstr)
        #cstr.properties.ndotin = 1e-3
        soln = solveIsothermal(sol.xin[8, j, i], cstr, soln.u)
        kk = LinRange(2, 10, 20)
        for k in kk
            cstr.properties.ndotin = 1 / 10^k
            soln = solveIsothermal(sol.xin[8, j, i], cstr, soln.u)
        end

        yield_fbr(soln.u, cstr)
        Y = zeros(40)


        for i = 1:40
            cstr.properties.yin = sol.xin[1:6, j, i]
            cstr.properties.Tc = sol.xin[8, j, i]




            #sol=solveIsothermal(200, cstr)#, x0=vcat(cstr.properties.yin, 0.5, cstr.properties.Tc, cstr.properties.ndotin))

            soln = solveIsothermal(sol.xin[8, j, i], cstr, soln.u)

            Y[i] = yield_fbr(soln.u, cstr)
        end
        Ym[:, j] = Y
    end
    return Ym
end

function calc_improve(solSS1, sol1)

    splf = Spline1D(sol1.obj2, sol1.obj1, k=1)
    sply = Spline1D(sol1.obj1[end:-1:1], sol1.obj2[end:-1:1], k=1)
    splssf = Spline1D(solSS1.obj2, solSS1.obj1, k=1)
    splssy = Spline1D(solSS1.obj1[end:-1:1], solSS1.obj2[end:-1:1], k=1)


    return [], [], splf, splssf, sply, splssy
end


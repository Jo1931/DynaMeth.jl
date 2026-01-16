function calcRHSbif(x, xin, mcat, A, reac::Reactor)
    # unpack parameter
    for i in eachindex(x[1:end-1])
        if x[i] < 0
            x[i] = 0
        end
    end
    @unpack kinetic, properties = reac
    @unpack cp_mix, k, H, P, cp_cat, Tc = properties
    @unpack parameter, rates, stoich, catalyst = kinetic


    y = x[1:7]
    T = x[8]
    ndotout = x[9]

    yin = xin[1:6]
    Tin = xin[8]
    ndotin = xin[9]
    gamma = reac.kinetic.prefaktor(y[7])
    r = rates(y, T, P / 1e5, parameter.beta, parameter) .* gamma



    RHS16 = ndotin * yin[1:6] - y[1:6] * ndotout + (mcat * stoich * r)
    RHS7 = catalyst(y, T, parameter.k, parameter.gibbs, parameter)
    RHS8 = ndotin * cp_mix * (Tin - T) - k * A * (T - Tc) - mcat * sum(H[i] * r[i] for i = 1:3)
    RHS9 = -ndotout + ndotin + mcat * (sum(stoich * r))

    RHS = [RHS16;  # component balances
        RHS7; #catalyst
        RHS8;#temperarture
        RHS9] #flow rate];

    return RHS
end
### Single CSTR with reflux
function calcRHSsingle(x, reac::Reactor)
    # unpack parameter
    for i in eachindex(x)
        if x[i] < 0
            x[i] = 0
        end
    end
    @unpack kinetic, properties = reac
    @unpack ndotin, cp_mix, k, A, mcat, H, yin, P, cp_cat, Vgas, reflux, Trec, Tc = properties
    @unpack parameter, rates, stoich, catalyst = kinetic


    y = x
    T = y[8]
    if Trec == "Tc"
        Trec = Tc
    else
        Trec = 500
    end
    xin = vcat(yin, 0, Tc, ndotin)
    RHSsingle = calcRHSbif(x, xin, mcat, A, reac::Reactor)

    yup = [0, 1, 1, 1, 0, 1] .* y[1:6]
    ydown = [1, 0, 0, 0, 1, 0] .* y[1:6]
    #RHS
    nup = y[9] * yup
    ndown = y[9] * ydown
    nrec = nup * reflux
    nwaste = nup * (1 - reflux)


    RHS16 = RHSsingle[1:6] + nrec
    RHS7 = RHSsingle[7]
    RHS8 = RHSsingle[8] + sum(nrec[i] for i = 1:6) * cp_mix * (Trec - T)
    RHS9 = RHSsingle[9] + sum(nrec[i] for i = 1:6)

    RHS = [RHS16;  # component balances
        RHS7; #catalyst
        RHS8;#temperarture
        RHS9] #flow rate];

    return RHS
end
### Cascade CSTR with reflux ###
function calcRHScascade(x, reac::Reactor)
    # unpack parameter
    for i in eachindex(x[1:end-1])
        if x[i] < 0
            x[i] = 0
        end
    end
    @unpack kinetic, properties = reac
    @unpack ndotin, cp_mix, k, A, mcat, H, yin, P, cp_cat, Vgas, reflux, Trec, Tc, num_casc = properties
    @unpack parameter, rates, stoich, catalyst = kinetic

    #xm = reshape(x, num_casc, 9)
    if Trec == "Tc"
        Trec0 = Tc
    elseif Trec == "Tend"
        Trec0 = x[8+(num_casc-1)*9]
    else
        Trec0 = 500
    end
    yup = [0, 1, 1, 1, 0, 1] .* x[(1:6).+(num_casc-1)*9]
    #RHS
    nup = x[9+(num_casc-1)*9] * yup
    nrec = nup * reflux
    nout = x[9+(num_casc-1)*9]
    Tin0 = Tc
    nges = (sum(nrec) + ndotin + nout)
    #Tin = (Tin0 * ndotin + Tin0 * sum(nrec)) / nges + Trec0 * sum(nout) / nges
    epsi = 0.5
    Tin = Tin0 + epsi * (Trec0 - Tin0)
    # Tin = Tin0 * ndotin / nges + Trec0 * (nout + sum(nrec)) / nges
    xin = vcat(yin, 0, Tin, ndotin)
    RHSsingle = calcRHSbif(x[1:9], xin, mcat / num_casc, A / num_casc, reac::Reactor)





    RHS16 = RHSsingle[1:6] + nrec
    RHS7 = RHSsingle[7]
    RHS8 = RHSsingle[8] + sum(nrec[i] for i = 1:6) * cp_mix * (Tin - x[8])
    RHS9 = RHSsingle[9] + sum(nrec[i] for i = 1:6)

    RHS = [RHS16;  # component balances
        RHS7;      #catalyst
        RHS8;      #temperarture
        RHS9]      #flow rate];
    if num_casc > 1
        for i = 2:num_casc
            RHS = vcat(RHS, calcRHSbif(x[(1:9).+(i-1)*9], x[(1:9).+(i-2)*9], mcat / num_casc, A / num_casc, reac::Reactor))
        end
    end
    return RHS
end
function calcRHSCascadeIsothermal(x, reac::Reactor)
    RHS = calcRHScascade(x, reac::Reactor)
    for i = 1:reac.properties.num_casc
        RHS[8+((i-1)*9)] = x[7*reac.properties.num_casc+i] - reac.properties.Tc
    end
    return RHS
end
function calcRHSbifIsothermal(x, reac::Reactor)
    RHS = calcRHSsingle(x, reac::Reactor)
    RHS[8] = x[8] - reac.properties.Tc
    return RHS
end
function calcRHSsingleArc(x, reac::Reactor)
    @unpack Tcprev, xprev, Tcpred, xpred, xpreprev, Tcpreprev, ds = reac.methods.predictor
    reac.properties.Tc = x[10]
    RHS = calcRHSsingle(x, reac::Reactor)
    Tc = x[10]
    arc = sum((x[i] - xprev[i])^2 for i = 1:9) + (Tc - Tcprev)^2 - ds^2
    return vcat(RHS, arc)
end
function calcRHSCascadeArc(x, reac::Reactor)
    @unpack Tcprev, xprev, Tcpred, xpred, xpreprev, Tcpreprev, ds = reac.methods.predictor
    reac.properties.Tc = x[end]
    RHS = calcRHScascade(x, reac::Reactor)
    Tc = x[end]
    arc = sum((x[i] - xprev[i])^2 for i = 1:length(x)-1) + (Tc - Tcprev)^2 - ds^2
    return vcat(RHS, arc)
end





function calc_ndotout(x, reac::Reactor)
    # unpack parameter
    for i in eachindex(x)
        if x[i] < 0
            x[i] = 0
        end
    end
    @unpack kinetic, properties = reac
    @unpack cp_mix, k, H, yin, P, cp_cat, reflux, Trec, Tc, ndotin, mcat, A = properties
    @unpack parameter, rates, stoich, catalyst = kinetic


    y = x[1:7]
    T = x[8]

    ndotin
    gamma = reac.kinetic.prefaktor(y[7])
    r = rates(y, T, P / 1e5, parameter.beta, parameter) .* gamma

    RHS = ndotin + mcat * (sum(stoich * r))



    return RHS
end
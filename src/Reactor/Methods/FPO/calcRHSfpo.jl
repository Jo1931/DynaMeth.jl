function calcRHSfpo(x, reac::Reactor)
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
    end

    Tin = Tc

    r = rates(y, T, P * 1e-5, parameter)

    #LHS
    yup = [0, 1, 1, 1, 0, 1] .* y[1:6]
    ydown = [1, 0, 0, 0, 1, 0] .* y[1:6]
    #RHS
    nup = y[9] * yup
    ndown = y[9] * ydown
    nrec = nup * reflux
    nwaste = nup * (1 - reflux)


    RHS16 = ndotin * yin[1:6] + nrec - y[1:6] * y[9] + (mcat * stoich * r)
    RHS7 = catalyst(y, T, parameter)
    RHS8 = ndotin * cp_mix * (Tin - T) + sum(nrec[i] for i = 1:6) * cp_mix * (Trec - T) - k * A * (T - Tc) - mcat * sum(H[i] * r[i] for i = 1:3)
    RHS9 = -y[9] + ndotin + sum(nrec[i] for i = 1:6) + mcat * (sum(stoich * r))

    RHS = [RHS16;  # component balances
        RHS7; #catalyst
        RHS8;#temperarture
        RHS9] #flow rate];

    return RHS
end
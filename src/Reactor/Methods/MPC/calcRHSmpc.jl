function calcRHSmpc(y, p, t0=0; T=p.properties.T, Pe=p.properties.P * 1e-5, yin=p.properties.yin, nin=p.properties.ndotin, activity=p.properties.activity)
    t = t0 - p.methods.deadtime
    mcat = p.properties.mcat
    @unpack stoich, prefaktor = p.kinetic
    @unpack beta, gibbs, k = p.kinetic.parameter
    gamma = prefaktor(y[7])
    r = activity * p.kinetic.rates(y, T, Pe, beta, p.kinetic.parameter) .* gamma .* p.methods.para
    #RHS

    ninv = calc_nin(p.methods.disturb, p.methods.input, yin, nin, t)


    #noutv = y[1:6] * y[8]
    ninn = sum(ninv)

    RHS16 = ninv - y[8] * y[1:6] + mcat * stoich * r
    RHS7 = p.kinetic.catalyst(y, T, k, gibbs, p.kinetic.parameter)
    RHS8 = -y[8] + ninn + mcat * (sum(stoich * r))
    RHS = [RHS16; RHS7; RHS8]
    return RHS
end

function calc_nin(dist, inlet, yin, nin, t)
    ninv = [nin * yin[1];
        inlet.co2(t);
        inlet.co(t);
        nin * yin[4] * (1 .+ dist(t));
        nin * yin[5:6]]
    return ninv
end

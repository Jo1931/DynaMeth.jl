function calc_x0(tv, cstr_plant)
    dy0 = zero(y0)   # Startwerte für die Ableitungen (häufig Null)
    differential_vars = [true, true, true, true, true, true, true, false]
    tend = tv[end]
    tspan = (0.0, tend)     # Zeitintervall
    prob = DAEProblem(dae_system_mpc, dy0, y0, tspan, cstr_plant, differential_vars=differential_vars, saveat=tv)
    sol = solve(prob, DFBDF())
    return sol[:, :]
end

function setInputCons(model, u, mean_N2_in, u0, cstr, ::Val{:CO_N2_NIN})
    N2 = cstr.methods.N2

    if !cstr.methods.optimize_temperature
        fix(u[10], cstr.properties.Tin; force=true)
        fix(u[11], cstr.properties.Tc; force=true)
    else
        @constraint(model, u[10] == u[11])
    end

    ts = LinRange(0, 2 * pi, 10000001)
    ratio = sum(cstr.methods.signal.(ts)) / length(ts)

    #@constraint(model, u[12] == cstr.properties.ndotin )#/ (1 + ratio * 0) * 100)
    fix(u[12], cstr.properties.ndotin; force=true)
    fix(u[13], 0; force=true)
    fix(u[14], 0; force=true)
    fix(u[15], 0; force=true)


    if cstr.methods.signal == cos
        @NLconstraint(model, c22, -u[9] * (u[6] * u[4] * cos(u[8] * pi) - 2) / 2 == N2)
        #@NLconstraint(model, abs(u[12]*p["N2"] - mean_N2_in) == 0)
        @constraint(model, u[5] * u[2] - u[6] * u[9] - u[14] * u[3] - u[13] * u[1] == 0)

        #@constraint(model,c22, u[9] == p["N2"])
    elseif cstr.methods.signal(0.1) == 0
        @constraint(model, c22, u[9] == N2)
        fix(u[4], 0; force=true)
        fix(u[5], 0; force=true)
        fix(u[6], 0; force=true)
        fix(u[7], 0.005; force=true)
        fix(u[8], 0; force=true)
    else
        @NLconstraint(model, abs(cstr.properties.ndotin * N2 - mean_N2_in) <= 1e-8)
        @constraint(model, u[5] * u[2] - u[6] * u[9] - u[14] * u[3] - u[13] * u[1] == 0)
    end
    @constraint(model, u[1] + u[2] >= 0.01)
    @constraint(model, u[1] + u[2] + u[3] + u[9] == 1)



    # @NLconstraint(model, n_in3+n_in2 >=0.01*u[12])

end
function setInputCons(model, u, mean_N2_in, u0, cstr, ::Val{:CO2_N2_NIN})
    N2 = cstr.methods.N2

    if !cstr.methods.optimize_temperature
        fix(u[10], cstr.properties.Tin; force=true)
        fix(u[11], cstr.properties.Tc; force=true)
    else
        @constraint(model, u[10] == u[11])
        #fix(u[10], cstr.properties.Tin; force=true)
    end

    ts = LinRange(0, 2 * pi, 10000001)
    ratio = sum(cstr.methods.signal.(ts)) / length(ts)

    #@constraint(model, u[12] == cstr.properties.ndotin )#/ (1 + ratio * 0) * 100)
    fix(u[12], cstr.properties.ndotin; force=true)
    fix(u[5], 0; force=true)
    fix(u[14], 0; force=true)
    fix(u[15], 0; force=true)
    fix(u[2], 0; force=true)


    if cstr.methods.signal == cos
        @NLconstraint(model, c22, -u[9] * (u[6] * u[4] * cos(u[8] * pi) - 2) / 2 == N2)
        #@NLconstraint(model, abs(u[12]*p["N2"] - mean_N2_in) == 0)
        @constraint(model, u[5] * u[2] - u[6] * u[9] - u[14] * u[3] - u[13] * u[1] == 0)

        #@constraint(model,c22, u[9] == p["N2"])
    elseif cstr.methods.signal(0.1) == 0
        @constraint(model, c22, u[9] == N2)
        fix(u[4], 0; force=true)
        fix(u[13], 0; force=true)
        fix(u[6], 0; force=true)
        fix(u[7], 0.005; force=true)
        fix(u[8], 0; force=true)
    else
        @NLconstraint(model, abs(cstr.properties.ndotin * N2 - mean_N2_in) <= 1e-5)
        @constraint(model, u[5] * u[2] - u[6] * u[9] - u[14] * u[3] - u[13] * u[1] == 0)
    end
    @constraint(model, u[1] + u[2] >= 0.01)
    @constraint(model, u[1] + u[2] + u[3] + u[9] == 1)



    # @NLconstraint(model, n_in3+n_in2 >=0.01*u[12])

end
function setInputCons(model, u, mean_N2_in, u0, cstr, ::Val{:METHANE})
    N2 = cstr.methods.N2

    if !cstr.methods.optimize_temperature
        fix(u[10], cstr.properties.Tin; force=true)
        fix(u[11], cstr.properties.Tc; force=true)
    else
        #@constraint(model, u[10] == u[11])
        fix(u[10], cstr.properties.Tin; force=true)
    end

    fix(u[2], 0; force=true)
    ts = LinRange(0, 2 * pi, 10000001)
    ratio = sum(cstr.methods.signal.(ts)) / length(ts)

    #@constraint(model, u[12] == cstr.properties.ndotin )#/ (1 + ratio * 0) * 100)
    fix(u[12], cstr.properties.ndotin; force=true)
    fix(u[5], 0; force=true)
    #fix(u[14], 0; force=true)
    fix(u[15], 0; force=true)
    fix(u[2], 0; force=true)

    #@constraint(model, u[3] == (1 - N2) / 5 * 4)
    #@constraint(model, u[3] == 4 * u[1])

    #fix(u[6], 1.0; force=true)
    if cstr.methods.signal == cos
        @NLconstraint(model, c22, -u[9] * (u[6] * u[4] * cos(u[8] * pi) - 2) / 2 == N2)
        #@NLconstraint(model, -u[3] * (u[14] * u[4] * cos(u[8] * pi) - 2) / 2 == 0.68)
        #@NLconstraint(model, c22, -u[3] * (u[6] * u[4] * cos(u[8] * pi) - 2) / 2 == N2)
        #@NLconstraint(model, abs(u[12]*p["N2"] - mean_N2_in) == 0)
        @constraint(model, u[5] * u[2] - u[6] * u[9] - u[14] * u[3] - u[13] * u[1] == 0)
        #@NLconstraint(model, abs(cstr.properties.ndotin * N2 - mean_N2_in) <= 1e-5)
        #@NLconstraint(model, abs(cstr.properties.ndotin * 0.68 - model[:mean_H2_in]) <= 1e-5)
        #@constraint(model,c22, u[9] == p["N2"])
    elseif cstr.methods.signal(0.1) == 0
        @constraint(model, c22, u[9] == N2)
        fix(u[4], 0; force=true)
        fix(u[13], 0; force=true)
        fix(u[14], 0; force=true)
        fix(u[6], 0; force=true)
        fix(u[7], 0.005; force=true)
        fix(u[8], 0; force=true)
        #@constraint(model, 4 * u[1] == u[3])
    else
        @NLconstraint(model, abs(cstr.properties.ndotin * N2 - mean_N2_in) <= 1e-5)
        @constraint(model, u[5] * u[2] - u[6] * u[9] - u[14] * u[3] - u[13] * u[1] == 0)
    end
    @constraint(model, u[1] + u[2] >= 0.01)
    @constraint(model, u[1] + u[2] + u[3] + u[9] == 1)



    # @NLconstraint(model, n_in3+n_in2 >=0.01*u[12])

end
function setInputCons(model, u, mean_N2_in, u0, cstr, ::Val{:NONE})

    #@constraint(model, u[10] == u[11])
    fix(u[1], u0[1]; force=true)
    fix(u[2], u0[2]; force=true)
    fix(u[3], u0[3]; force=true)
    fix(u[4], u0[4]; force=true)
    fix(u[5], u0[5]; force=true)
    fix(u[6], u0[6]; force=true)
    fix(u[7], u0[7]; force=true)
    fix(u[8], u0[8]; force=true)
    fix(u[9], u0[9]; force=true)
    fix(u[10], cstr.properties.Tin; force=true)
    fix(u[11], cstr.properties.Tc; force=true)
    fix(u[12], cstr.properties.ndotin; force=true)
    fix(u[13], u0[13]; force=true)
    fix(u[14], u0[14]; force=true)
    fix(u[15], u0[15]; force=true)


end


function setInputCons(model, u, mean_N2_in, u0, cstr, ::Val{:CO2_H2_NIN_NOCO})
    N2 = cstr.methods.N2
    @constraint(model, u[10] == u[11])
    fix(u[11], cstr.properties.Tc; force=true)
    fix(u[12], cstr.properties.ndotin; force=true)
    fix(u[5], 0; force=true)
    fix(u[6], 0; force=true)
    fix(u[15], 0; force=true)
    fix(u[2], 0; force=true)
    @constraint(model, c22, u[9] == N2)


    if cstr.methods.signal == cos

        @constraint(model, u[5] * u[2] - u[6] * u[9] - u[14] * u[3] + u[13] * u[1] == 0)

    else
        cstr.methods.signal(0.1) == 0
        fix(u[4], 0; force=true)
        fix(u[13], 0; force=true)
        fix(u[14], 0; force=true)
        fix(u[7], 0.005; force=true)
        fix(u[8], 0; force=true)
    end
    @constraint(model, u[1] + u[2] >= 0.01)
    @constraint(model, u[1] + u[2] + u[3] + u[9] == 1)



    # @NLconstraint(model, n_in3+n_in2 >=0.01*u[12])

end

function setInputCons(model, u, mean_N2_in, u0, cstr, op::Symbol)
    return setInputCons(model, u, mean_N2_in, u0, cstr, Val(op))
end

function daeFuncMess(out, dy, y, p, t)
    reac = p[2]
    spl = p[1]

    ## calc inlets

    xin1 = 0
    xin2 = spl[1](t)
    xin3 = spl[2](t)
    xin4 = spl[3](t)
    xin5 = 0.0005
    xin6 = spl[4](t) - xin5
    xin7 = 0
    yin = [xin1; xin2; xin3; xin4; xin5; xin6; xin7]
    nin = spl[5](t)
    T = spl[6](t)
    P = spl[7](t)
    activity = reac.properties.activity
    param = get_kinetic_param(reac.kinetic.parameter)
    RHS = calcRHSexp(y, reac, T, P * 1e-5, yin, nin, param, activity)
    LHS = calcLHSexp(y, reac)
    out .= LHS * dy - RHS
end

function odeFuncMess(out, y, p, t)
    reac = p[2]
    spl = p[1]

    ## calc inlets

    xin1 = 0
    xin2 = spl[1](t)
    xin3 = spl[2](t)
    xin4 = spl[3](t)
    xin5 = 0
    xin6 = spl[4](t)
    xin7 = 0
    reac.properties.yin = [xin1; xin2; xin3; xin4; xin5; xin6; xin7]
    reac.properties.ndotin = spl[5](t)
    reac.properties.T = spl[6](t)
    reac.properties.P = spl[7](t)

    x = y
    for i = 1:7
        if x[i] <= 0
            x[i] = 0
        end
    end
    RHS = calcRHSexp(x, reac)
    LHS = calcLHSexp(x, reac)
    out .= inv(LHS) * RHS
end



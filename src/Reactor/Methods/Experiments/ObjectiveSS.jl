struct ObjectiveSS end
function (res::ObjectiveSS)(x, cstr::T=cstr) where {T<:Union{Reactor,ComponentVector}}
    #x = unwrappedParam(u, cstr)
    #setDecisionVar!(cstr, x)

    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    u = unwrappedParam(x, cstr)
    xsim = solveSS(cstr, i, param=u)
    xmess = cstr.methods.xout_mess[i]
    out = sum(((xsim[1:3] - xmess[1:3]) ./ (xmess[1:3] .+ 1e-1)) .^ 2)
    #out = sum(((xsim[1:3] - xmess[1:3]) .^ 2))
    for i in experimental_set[2:end]
        xsims = solveSS(cstr, i, param=u)
        xmesss = cstr.methods.xout_mess[i]
        out = out + sum(((xsims[1:3] - xmesss[1:3]) ./ (xmesss[1:3] .+ 1e-1)) .^ 2)
        #out = out + sum(((xsims[1:3] - xmesss[1:3]) .^ 2))
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end
    #penalty = penaltyFunc(cstr)
    if cstr.methods.idx == "Dyn"
        dyn = cstr.methods.objectiveDyn(x, cstr)
    else
        dyn = 0
    end
    out = out .+ dyn * 50

    cstr.methods.xout_sim = xsim
    res = out #+ regularization(u)
end
#function (res::ObjectiveSS)(cstr::Reactor)
#    experimental_set = cstr.methods.subset
#    i = experimental_set[1]
#    u = get_kinetic_param(cstr.kinetic.parameter)
#    xsim = solveSS(cstr, i, param=u)
#    xmess = cstr.methods.xout_mess[i]
#    out = sum(((xsim[1:3] - xmess[1:3]) ./ (xmess[1:3] .+ 1e-1)) .^ 2)
#out = sum(((xsim[1:3] - xmess[1:3]) .^ 2))
#    for i in experimental_set[2:end]
#        xsims = solveSS(cstr, i, param=u)
#        xmesss = cstr.methods.xout_mess[i]
#        out = out + sum(((xsims[1:3] - xmesss[1:3]) ./ (xmesss[1:3] .+ 1e-1)) .^ 2)
#        #out = out + sum(((xsims[1:3] - xmesss[1:3]) .^ 2))
#       xsim = hcat(xsim, xsims)
#        xmess = hcat(xmess, xmesss)
#    end
#penalty = penaltyFunc(cstr)
#    if cstr.methods.idx == "Dyn"
#       dyn = cstr.methods.objectiveDyn(x, cstr)
#   else
#       dyn = 0
#   end
#   out = out .+ dyn * 50
#   cstr.methods.xout_sim = xsim
#   res = out
#end
function penaltyFunc(cstr::Reactor, p::Seidel)
    penalty = minimum([0, cstr.kinetic.parameter.beta[2]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[4]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[6]]) .^ 2
    for i = 7:14
        penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[i]]) .^ 2
    end
    for i = 1:2
        penalty = penalty + minimum([0, cstr.kinetic.parameter.k[i]]) .^ 2
    end
    for i = 1:2
        penalty = penalty + minimum([0, cstr.kinetic.parameter.gibbs[i]]) .^ 2
    end
    penalty = penalty + minimum([0, 1 - cstr.kinetic.parameter.beta[10]]) .^ 2
    return penalty
end
function penaltyFunc(cstr::Reactor)
    penaltyFunc(cstr, cstr.kinetic.ident)
end
function penaltyFunc(cstr::Reactor, p::AutoCatSeidel)
    penalty = minimum([0, cstr.kinetic.parameter.beta[2]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[4]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[6]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[16]]) .^ 2
    for i = 7:14
        penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[i]]) .^ 2
    end
    for i = 1:2
        penalty = penalty + minimum([0, cstr.kinetic.parameter.k[i]]) .^ 2
    end
    for i = 1:2
        penalty = penalty + minimum([0, cstr.kinetic.parameter.gibbs[i]]) .^ 2
    end
    penalty = penalty + minimum([0, 1 - cstr.kinetic.parameter.beta[10]]) .^ 2
    return penalty
end
function penaltyFunc(cstr::Reactor, p::T) where {T<:Union{AutoCat,AutoCatLumped}}
    penalty = minimum([0, cstr.kinetic.parameter.beta[2]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[4]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[6]]) .^ 2
    penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[8]]) .^ 2
    for i = 9:28
        penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[i]]) .^ 2
    end
    # for i = 1:2
    #     penalty = penalty + minimum([0, cstr.kinetic.parameter.k[i]]) .^ 2
    # end
    # for i = 1:2
    #     penalty = penalty + minimum([0, cstr.kinetic.parameter.gibbs[i]]) .^ 2
    # end
    # penalty = penalty + minimum([0, 1 - cstr.kinetic.parameter.beta[10]]) .^ 2
    return penalty
end
function penaltyFunc(cstr::Reactor, p::Nestler)
    return 0
end
function penaltyFuncGrad(x)
    penalty = zeros(18)
    for i in [2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        penalty[i] = minimum([0, x[i]]) * 2
    end
    return penalty
end
function regularization(beta)
    fac = 1e-8
    fac1 = 1e-8
    #fac * ((beta[2] - 49.28588684344555)^2 + (beta[4] - 22.1373)^2 + (beta[6] - 11.1406)^2 + (beta[8] - 32.8078)^2) + fac1 * sum(beta[i]^2 for i = 9:28) #+ (beta[2] - 49.28588684344555)^2 + (beta[4] - 20.1679)^2 + (beta[6] - 13.3051)^2 + (beta[8] - 35.7427)^2)
    #fac * ((beta[2] - 49.28588684344555)^2 + 30 * (beta[4] - 23.602430354530163)^2 + ((beta[6] - 13.037502048475364))^2 + (beta[8] - 29.03433353812334)^2) + fac1 * sum(beta[i]^2 for i = 9:28) #+ 10 * beta[28]^2#+ (beta[2] - 49.28588684344555)^2 + (beta[4] - 20.1679)^2 + (beta[6] - 13.3051)^2 + (beta[8] - 35.7427)^2)
    #fac * ((beta[2] - 49.28588684344555)^2 + (beta[4] - 20.1679)^2            + (beta[6] - 13.3051)^2            + (beta[8] - 35.7427)^2) + fac1 * sum(beta[i]^2 for i = 9:28) #+ (beta[2] - 49.28588684344555)^2 + (beta[4] - 20.1679)^2 + (beta[6] - 13.3051)^2 + (beta[8] - 35.74

    beta0 = [-5.001, 26.455, -3.145, 1.5308, -4.4526, 15.615, 1.1064, 0.0, 0.0, 0.0, 0.14969, 0.0, 0.062881, 0.0, 79.174 * 1e-4, 0.188 * 1e-4, 0.335747716622, 2.184147818E+1]

    #k = [79.174, 0.188] * 1e-4
    #fac * sum((beta[i] - beta0[i])^2 for i = 1:6) + fac1 * sum((beta[i] - beta0[i])^2 for i = 7:18)

    #    fac * ((beta[2] - 48.45707829944559)^2 + 30 * (beta[4] - 31.60616568134138)^2 + ((beta[6] - 13.654562896921291))^2 + (beta[8] - 34.82038193189341)^2) + fac1 * sum(beta[i]^2 for i = 9:28) #+ 10 * beta[28]^2#+ (beta[2] - 49.28588684344555)^2 + (beta[4] - 20.1679)^2 + (beta[6] - 13.3051)^2 + (beta[8] - 35.7427)^2)
    # 0
    # 1e-6 * sum((beta[i])^2 for i = 7:14)
    #beta0 = [4.779 * 1e10, 1.142 * 1e5, 3.529 * 1e14, 1.448 * 1e5, 3.052, 9.695 * 1e2, 3.545 * 1e1]
    1e-20 * sum((beta[i] - beta0[i])^2 for i = 1:18)
end
#21.6349  6.63129  17.1557  17.323  12.0793

function obj_Rexp(x, cstr)
    u = unwrappedParam(x, cstr)
    experimental_set = cstr.methods.subset
    i = experimental_set[1]

    P = cstr.methods.xin_mess[i][9]
    T = cstr.methods.xin_mess[i][8]
    y = cstr.methods.xout_mess[i][1:6]
    nin = cstr.methods.xin_mess[i][7]
    yin = cstr.methods.xin_mess[i][1:6]
    Rsim = cstr.kinetic.stoich * cstr.kinetic.rates(y, T, P * 1e-5, u, cstr.kinetic.parameter)
    gamma = cstr.methods.xin_mess[i][6] / cstr.methods.xout_mess[i][6]
    Rexp = -(yin .- y * gamma) * nin / cstr.properties.mcat
    ob = sum((Rexp[1:3] .- Rsim[1:3]) .^ 2 ./ (Rexp[1:3] .+ 1e-1))
    for i in experimental_set[2:end]
        P = cstr.methods.xin_mess[i][9]
        T = cstr.methods.xin_mess[i][8]
        y = cstr.methods.xout_mess[i][1:6]
        nin = cstr.methods.xin_mess[i][7]
        yin = cstr.methods.xin_mess[i][1:6]
        Rsim = cstr.kinetic.stoich * cstr.kinetic.rates(y, T, P * 1e-5, u, cstr.kinetic.parameter)
        gamma = cstr.methods.xin_mess[i][6] / cstr.methods.xout_mess[i][6]
        Rexp = -(yin .- y * gamma) * nin / cstr.properties.mcat
        ob = ob + sum((Rexp[1:3] - Rsim[1:3]) .^ 2 / (Rexp[1:3] .+ 1e-1))
    end
    return ob
end



function ObjectiveSS(x, p)
    #x = unwrappedParam(u, cstr)
    #setDecisionVar!(cstr, x)

    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    u = unwrappedParam(x, cstr)
    xsim = solveSS(cstr, i, param=u)
    xmess = cstr.methods.xout_mess[i]
    out = sum(((xsim[1:3] - xmess[1:3]) ./ (xmess[1:3] .+ 1e-1)) .^ 2)
    #out = sum(((xsim[1:3] - xmess[1:3]) .^ 2))
    for i in experimental_set[2:end]
        xsims = solveSS(cstr, i, param=u)
        xmesss = cstr.methods.xout_mess[i]
        out = out + sum(((xsims[1:3] - xmesss[1:3]) ./ (xmesss[1:3] .+ 1e-1)) .^ 2)
        #out = out + sum(((xsims[1:3] - xmesss[1:3]) .^ 2))
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end
    #penalty = penaltyFunc(cstr)

    cstr.methods.xout_sim = xsim
    res = out + regularization(u)
end


function (res::ObjectiveSS)(cstr::Reactor)
    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    u = get_kinetic_param(cstr.kinetic.parameter)
    xsim = solveSS(cstr, i, param=u)
    xmess = cstr.methods.xout_mess[i]
    out = sum(((xsim[1:3] - xmess[1:3]) ./ (xmess[1:3] .+ 1e-1)) .^ 2)
    #out = sum(((xsim[1:3] - xmess[1:3]) .^ 2))
    for i in experimental_set[2:end]
        xsims = solveSS(cstr, i, param=u)
        xmesss = cstr.methods.xout_mess[i]
        out = out + sum(((xsims[1:3] - xmesss[1:3]) ./ (xmesss[1:3] .+ 1e-1)) .^ 2)
        #out = out + sum(((xsims[1:3] - xmesss[1:3]) .^ 2))
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end
    #penalty = penaltyFunc(cstr)
    if cstr.methods.idx == "Dyn"
        dyn = cstr.methods.objectiveDyn(x, cstr)
    else
        dyn = 0
    end
    out = out .+ dyn * 50
    cstr.methods.xout_sim = xsim
    res = out
end


function RMSE(cstr::Reactor)
    experimental_set = cstr.methods.subset
    i = experimental_set[1]
    u = get_kinetic_param(cstr.kinetic.parameter)
    xsim = solveSS(cstr, i, param=u)
    xmess = cstr.methods.xout_mess[i]
    out = sum(((xsim[1:3] - xmess[1:3])) .^ 2)
    #out = sum(((xsim[1:3] - xmess[1:3]) .^ 2))
    for i in experimental_set[2:end]
        xsims = solveSS(cstr, i, param=u)
        xmesss = cstr.methods.xout_mess[i]
        out = out + sum(((xsims[1:3] - xmesss[1:3])) .^ 2)
        #out = out + sum(((xsims[1:3] - xmesss[1:3]) .^ 2))
        xsim = hcat(xsim, xsims)
        xmess = hcat(xmess, xmesss)
    end
    #penalty = penaltyFunc(cstr)
    if cstr.methods.idx == "Dyn"
        dyn = cstr.methods.objectiveDyn(x, cstr)
    else
        dyn = 0
    end
    out = out .+ dyn * 50
    cstr.methods.xout_sim = xsim
    res = sqrt(out / length(experimental_set))
end
function getDecisionVar(cstr::Reactor, p::Seidel)
    return vcat(cstr.kinetic.parameter.beta, cstr.kinetic.parameter.k, cstr.kinetic.parameter.gibbs)
end
function getDecisionVar(cstr::Reactor, p::AutoCatSeidel)
    return vcat(cstr.kinetic.parameter.beta, cstr.kinetic.parameter.k, cstr.kinetic.parameter.gibbs)
end
function getDecisionVar(cstr::Reactor, p::AutoCat)
    return vcat(cstr.kinetic.parameter.beta, cstr.kinetic.parameter.k, cstr.kinetic.parameter.gibbs)
end
function getDecisionVar(cstr::Reactor, p::AutoCatLumped)
    return vcat(cstr.kinetic.parameter.beta, cstr.kinetic.parameter.k, cstr.kinetic.parameter.gibbs)
end
function getDecisionVar(cstr::Reactor, p::Nestler)
    return vcat(cstr.kinetic.parameter.beta)
end
function setDecisionVar!(cstr::Reactor, u, p::Seidel)
    cstr.kinetic.parameter.beta = u[1:14]
    cstr.kinetic.parameter.k = u[15:16]
    cstr.kinetic.parameter.gibbs = u[17:18]
    nothing
end
function setDecisionVar!(cstr::Reactor, u, p::AutoCatSeidel)
    cstr.kinetic.parameter.beta[1:16] = u[1:16]
    #cstr.kinetic.parameter.k = u[17:18]
    #cstr.kinetic.parameter.gibbs = u[19:20]
    nothing
end
function setDecisionVar!(cstr::Reactor, u, p::T) where {T<:Union{AutoCat,AutoCatLumped}}
    cstr.kinetic.parameter.beta[1:length(u)] = u
    # cstr.kinetic.parameter.k = u[29:30]
    # cstr.kinetic.parameter.gibbs = u[31:32]
    nothing
end
function setDecisionVar!(cstr::Reactor, u, p::Nestler)
    cstr.kinetic.parameter.beta = u[1:9]
    #cstr.kinetic.parameter.k = u[29:30]
    #cstr.kinetic.parameter.gibbs = u[31:32]
    nothing
end


function getDecisionVar(cstr::Reactor)
    getDecisionVar(cstr::Reactor, cstr.kinetic.ident)
end
function setDecisionVar!(cstr::Reactor, u)
    setDecisionVar!(cstr::Reactor, u, cstr.kinetic.ident)
end
function wrappedParam(u, p::Seidel)
    uw = copy(u)
    for i in [2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    return uw
end
function unwrappedParam(uw, p::Seidel)
    u = copy(uw)
    for i in [2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    return u
end
function wrappedParam(u, cstr::Reactor)
    wrappedParam(u, typeof(cstr.kinetic.rates))
end
function unwrappedParam(u, cstr::Reactor)
    unwrappedParam(u, typeof(cstr.kinetic.rates))
end
function wrappedParam(u, p::Nestler)
    uw = copy(u)
    typical_x = [5.441 * 1e-4, -45458, 24.701, -54970, 33.321 * 1e-18, 109959, 8.262 * 1e-6, 6.430 * 1e-14, 119570]
    return uw ./ typical_x #sign.(uw) .* abs.(uw) .^ (1 / 3)
end
function unwrappedParam(uw, p::Nestler)
    u = copy(uw)
    typical_x = [5.441 * 1e-4, -45458, 24.701, -54970, 33.321 * 1e-18, 109959, 8.262 * 1e-6, 6.430 * 1e-14, 119570]
    return u .* typical_x#u .^ (3)
end
function wrappedParam(u, ::T) where {T<:Union{Type{ReactionRatesAutoCat},Type{ReactionRatesAutoCatR2},Type{ReactionRatesAutoCatR4},Type{ReactionRatesAutoCatWithoutCO},Type{ReactionRatesAutoCatWithoutCOR2}}}
    uw = copy(u)
    for i in [2, 4, 6, 8]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    for i in 9:length(u)#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    return uw
end
function unwrappedParam(uw, ::T) where {T<:Union{Type{ReactionRatesAutoCat},Type{ReactionRatesAutoCatR2},Type{ReactionRatesAutoCatR4},Type{ReactionRatesAutoCatWithoutCO},Type{ReactionRatesAutoCatWithoutCOR2}}}
    u = copy(uw)
    for i in [2, 4, 6, 8]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    for i in 9:length(u)#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    return u
end
function wrappedParam(u, p::AutoCatLumped)
    uw = copy(u)
    for i in [2, 4, 6, 8]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    for i in 9:length(u)#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    return uw
end
function unwrappedParam(uw, p::AutoCatLumped)
    u = copy(uw)
    for i in [2, 4, 6, 8]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    for i in 9:length(u)#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    return u
end
function wrappedParam(u, ::Type{ReactionRatesSeidel})
    uw = copy(u)
    for i in [2, 4, 6]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    for i in 7:18#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    return uw
end
function unwrappedParam(uw, ::Type{ReactionRatesSeidel})
    u = copy(uw)
    for i in [2, 4, 6]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    for i in 7:18#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    return u
end
function wrappedParam(u, ::Type{ReactionRatesNestler})
    uw = copy(u)
    #beta = [5.441 * 1e-4, -45458, 24.701, -54970, 33.321 * 1e-18, 109959, 8.262 * 1e-6, 6.430 * 1e-14, 119570, 1, 1, 1, 1]
    beta = [-17.9733, 10.4514, -9.4315, 12.6383, -12.6593, -25.2810, 8.262 * 1e-6, -2.8845, -27.4907, 1, 1, 1, 1]

    return uw ./ beta
end
function unwrappedParam(uw, ::Type{ReactionRatesNestler})
    u = copy(uw)
    #beta = [5.441 * 1e-4, -45458, 24.701, -54970, 33.321 * 1e-18, 109959, 8.262 * 1e-6, 6.430 * 1e-14, 119570, 1, 1, 1, 1]
    beta = [-17.9733, 10.4514, -9.4315, 12.6383, -12.6593, -25.2810, 8.262 * 1e-6, -2.8845, -27.4907, 1, 1, 1, 1]

    return u .* beta
end
function wrappedParam(u, ::Type{ReactionRatesBussche})
    uw = copy(u)
    #beta = [1.07 * 1e-10, 36696, 1.22 * 1e5, -94765, 3.45338, 1.578 * 1e-3, 17197, 6.621e-16, 124119, 1, 1, 1, 1]
    beta = [-14.5213, -8.4369, -10.0759, 21.7877, 3.45338, -2.4978, -3.9538, -6.4147, -28.5366, 1, 1, 1, 1]
    return uw ./ beta
end
function unwrappedParam(uw, ::Type{ReactionRatesBussche})
    u = copy(uw)
    #beta = [1.07 * 1e-10, 36696, 1.22 * 1e5, -94765, 3.45338, 1.578 * 1e-3, 17197, 6.621e-16, 124119, 1, 1, 1, 1]
    beta = [-14.5213, -8.4369, -10.0759, 21.7877, 3.45338, -2.4978, -3.9538, -6.4147, -28.5366, 1, 1, 1, 1]
    return u .* beta
end
function wrappedParam(u, ::Type{ReactionRatesGraaf})
    uw = copy(u)
    beta = [4.89 * 1e7, -113000, 1.09 * 1e5, -87500, 9.64 * 1e11, -152900, 2.16 * 1e-5, 46800, 7.05 * 1e-7, 61700, 6.37 * 1e-9, 84000, 1, 1, 1, 1]
    return uw ./ beta
end
function unwrappedParam(uw, ::Type{ReactionRatesGraaf})
    u = copy(uw)
    beta = [4.89 * 1e7, -113000, 1.09 * 1e5, -87500, 9.64 * 1e11, -152900, 2.16 * 1e-5, 46800, 7.05 * 1e-7, 61700, 6.37 * 1e-9, 84000, 1, 1, 1, 1]
    return u .* beta
end
function wrappedParam(u, ::Type{ReactionRatesBrilman})
    uw = copy(u)
    for i in [2,4,5,6,7]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        uw[i] = log(u[i])
    end
    return uw
end
function unwrappedParam(uw, ::Type{ReactionRatesBrilman})
    u = copy(uw)
    for i in [2,4,5,6,7]#, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        u[i] = exp(uw[i])
    end
    return u
end

function set_kinetic_param!(cstr::Reactor, param)
    cstr.kinetic.parameter.beta = param[1:end-4]
    cstr.kinetic.parameter.k = param[end-3:end-2]
    cstr.kinetic.parameter.gibbs = param[end-1:end]
end
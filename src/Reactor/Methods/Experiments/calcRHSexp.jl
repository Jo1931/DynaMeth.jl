function calcRHSexp(y, p, T=p.properties.T, Pe=p.properties.P * 1e-5, yin=p.properties.yin, nin=p.properties.ndotin, param=get_kinetic_param(p.kinetic.parameter), activity=p.properties.activity)
    mcat = p.properties.mcat
    stoich = p.kinetic.stoich
    beta, k, gibbs = get_decision_variable(param)
    gamma = p.kinetic.prefaktor(y[7])
    r = activity * p.kinetic.rates(y, T, Pe, beta, p.kinetic.parameter) .* gamma
    #RHS
    RHS = [nin * (yin[1:6] - y[1:6]) + (mcat * (Matrix(I, 6, 6) - y[1:6] * ones(1, 6)) * stoich * r);
        p.kinetic.catalyst(y, T, k, gibbs, p.kinetic.parameter)]
    return RHS
end
function get_decision_variable(p::T) where {T<:Union{ComponentVector,NamedTuple}}
    @unpack beta, k, gibbs = p
    return beta, k, gibbs
end
function get_decision_variable(param::T) where {T<:Vector}
    beta = param[1:end-4]
    k = param[end-3:end-2]
    gibbs = param[end-1:end]
    return beta, k, gibbs
end


function calcHeatProduction(x, T, mcat, p::Reactor)
    @unpack properties, kinetic = p
    @unpack H, P = properties
    @unpack parameter, rates = kinetic
    gamma = p.kinetic.prefaktor(x[7])
    return -mcat * H' * (rates(x, T, P / 1e5, parameter.beta, parameter) .* gamma)
end
function calcHeatRemoval(T, Tin, Tc, ndotin, A, p::Reactor)
    @unpack properties = p
    @unpack cp_mix, k = properties
    return (ndotin * cp_mix) * (T - Tin) + k * A * (T - Tc)
end
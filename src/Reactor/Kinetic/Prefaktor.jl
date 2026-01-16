struct PrefaktorSeidel <: AbstractKineticComponent end
"""
    (J::PrefaktorSeidel)(phi)

Impact of catalyst on reaction rates as in Seidel et al. 2018
"""
function (J::PrefaktorSeidel)(phi)
    return [(1 - phi); (phi)^(2); phi / (1 - phi)]
end

struct PrefaktorNone <: AbstractKineticComponent end
"""
    (J::PrefaktorNone)(phi)

Impact of catalyst on reaction rates is not existing
"""
function (J::PrefaktorNone)(phi)
    return [1; 1; 1]
end
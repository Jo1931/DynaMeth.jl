"""
    Kinetic(stoich,parameter,rates,catalyst,jacobianLumped,prefaktor)

Parameter struct for all things related to the used kinetic for the reactor
"""
struct Kinetic <: AbstractReactor
    stoich
    parameter
    rates
    catalyst
    jacobianLumped
    prefaktor
end
"""
    Kinetic(p::AbstractKineticSetup)

Parameter struct constructor. p denotes setup for chosen kinetic. possible kinetics:

    Seidel()
    AutoCatSeidel()
    AutoCat()
    AutoCatLumped
    Nestler()
    Graaf()
    Bussche()
    Slotboom()
    Campos()
    Koschany()  
"""
function Kinetic(p::T) where {T<:AbstractKineticSetup}
    Kinetic(p.stoich, p.parameter, p.rates, p.catalyst, p.jacobian, p.prefaktor)
end
"""
    Seidel()

Setup for Methanol Kinetic: Seidel et al. 2018
    
Possible inputs:

    parameter       (default: KineticParameterSeidel())
    rates           (default: ReactionRatesSeidel())
    catalyst        (default: CatalystStateDynamic())
    jacobian        (default: JacobianLumpedSeidel())
    prefaktor       (default: PrefaktorSeidel())
"""
@kwdef struct Seidel <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterSeidel()
    rates = ReactionRatesSeidel()
    catalyst = CatalystStateDynamic()
    jacobian = JacobianLumpedSeidel()
    prefaktor = PrefaktorSeidel()
end
"""
    Seidel()

Setup for Methanol Kinetic: Seidel et al. 2018
    
Possible inputs:

    parameter       (default: KineticParameterSeidel())
    rates           (default: ReactionRatesSeidel())
    catalyst        (default: CatalystStateDynamic())
    jacobian        (default: JacobianLumpedSeidel())
    prefaktor       (default: PrefaktorSeidel())
"""
@kwdef struct Hybrid <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterHybrid()
    rates = ReactionRatesHybrid()
    catalyst = CatalystStateDynamic()
    jacobian = JacobianLumpedSeidel()
    prefaktor = PrefaktorNone()
end
"""
    AutoCatSeidel()

Setup for Methanol Kinetic: Seidel et al. 2018 with extention of autocatalytic step for co2 hydrogenation
    
Possible inputs:

    parameter       (default: KineticParameterAutoSeidel())
    rates           (default: ReactionRatesAutoSeidel())
    catalyst        (default: CatalystStateDynamic())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorSeidel())
"""
@kwdef struct AutoCatSeidel <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterAutoCatSeidel()
    rates = ReactionRatesAutoCatSeidel()
    catalyst = CatalystStateDynamic()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorSeidel()
end
"""
    AutoCat()

Setup for Methanol Kinetic: New developed by Kortuz with auto catalytic step for co2 hydrogenation
    
Possible inputs:

    parameter       (default: KineticParameterAutoCat())
    rates           (default: ReactionRatesAutoCat())
    catalyst        (default: CatalystStateDynamic())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct AutoCat <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 1; -1 0 -1; -2 -3 1; 0 1 -1; 0 0 0]
    parameter = KineticParameterAutoCat()
    rates = ReactionRatesAutoCat()
    catalyst = CatalystStateDynamic()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    AutoCatLumped()

Setup for lumped Methanol Kinetic: New developed by Kortuz with auto catalytic step for co2 hydrogenation
    
Possible inputs:

    parameter       (default: KineticParameterAutoLumped())
    rates           (default: ReactionRatesAutoLumped())
    catalyst        (default: CatalystStateDynamic())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct AutoCatLumped <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 1; -1 0 -1; -2 -3 1; 0 1 -1; 0 0 0]
    parameter = KineticParameterAutoCatLumped()
    rates = ReactionRatesAutoCatLumped()
    catalyst = CatalystStateDynamic()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    Nestler()

Setup for Methanol Kinetic: Nestler et al. 2020 (Graaf without CO-hydrogenation)
    
Possible inputs:

    parameter       (default: KineticParameterNestler())
    rates           (default: ReactionRatesNestler())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct Nestler <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterNestler()
    rates = ReactionRatesNestler()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    Graaf()

Setup for Methanol Kinetic: Graaf et al. 1988 
    
Possible inputs:

    parameter       (default: KineticParameterGraaf())
    rates           (default: ReactionRatesGraaf())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct Graaf <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterGraaf()
    rates = ReactionRatesGraaf()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    GraafRefit()

Setup for Methanol Kinetic: Graaf et al. 1988 
    
Possible inputs:

    parameter       (default: KineticParameterGraaf())
    rates           (default: ReactionRatesGraaf())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct GraafRefit <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterGraafRefit()
    rates = ReactionRatesGraafRefit()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    Bussche()

Setup for Methanol Kinetic: Vanden Bussche et al. 1996
    
Possible inputs:

    parameter       (default: KineticParameterBussche())
    rates           (default: ReactionRatesBussche())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct Bussche <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterBussche()
    rates = ReactionRatesBussche()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    Slotboom()

Setup for Methanol Kinetic: Slotboom et al. 2020
    
Possible inputs:

    parameter       (default: KineticParameterSlotboom())
    rates           (default: ReactionRatesSlottboom())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct Slotboom <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterSlotboom()
    rates = ReactionRatesSlotboom()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    Campos()

Setup for Methanol Kinetic: Campos et al. 2022
    
Possible inputs:

    parameter       (default: KineticParameterCampos())
    rates           (default: ReactionRatesCampos())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct Campos <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 1; -1 0 -1; -2 -3 1; 0 1 -1; 0 0 0]
    parameter = KineticParameterCampos()
    rates = ReactionRatesCampos()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    Brilman()

Setup for Methanol Kinetic: Brilman et al. 2025
    
Possible inputs:

    parameter       (default: KineticParameterCampos())
    rates           (default: ReactionRatesCampos())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct Brilman <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterBrilman()
    rates = ReactionRatesBrilman()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    BrilmanFull()

Setup for Methanol Kinetic: Brilman et al. 2025
    
Possible inputs:

    parameter       (default: KineticParameterCampos())
    rates           (default: ReactionRatesCampos())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct BrilmanFull <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterBrilmanFull()
    rates = ReactionRatesBrilmanFull()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
    Koschany()

Setup for Methane Kinetic: Koschany et al. 2016
    
Possible inputs:

    parameter       (default: KineticParameterKoschany())
    rates           (default: ReactionRatesKoschany())
    catalyst        (default: CatalystStateNone())
    jacobian        (default: JacobianLumpeNone())
    prefaktor       (default: PrefaktorNone())
"""
@kwdef struct Koschany <: AbstractKineticSetup
    stoich = [1 0 0; -1 0 0; 0 0 0; -4 0 0; 2 0 0; 0 0 0]
    parameter = KineticParameterKoschany()
    rates = ReactionRatesKoschany()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
"""
SeidelJump()

Setup for Methanol Kinetic: Seidel et al. 2018
    
Possible inputs:

    parameter       (default: KineticParameterSeidel())
    rates           (default: ReactionRatesSeidel())
    catalyst        (default: CatalystStateDynamic())
    jacobian        (default: JacobianLumpedSeidel())
    prefaktor       (default: PrefaktorSeidel())
"""
@kwdef struct SeidelJump <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterSeidel()
    rates = ReactionRatesSeidelJump()
    catalyst = CatalystStateDynamicJump()
    jacobian = JacobianLumpedSeidelJump()
    prefaktor = []
end
@kwdef struct KoschanyJump <: AbstractKineticSetup
    stoich = [1 0 0; -1 0 0; 0 0 0; -4 0 0; 2 0 0; 0 0 0]
    parameter = KineticParameterKoschany()
    rates = ReactionRatesKoschanyJump()
    catalyst = CatalystSteadyStateJump()
    jacobian = JacobianLumpedNoneJump()
    prefaktor = []
end
@kwdef struct HybridJump <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterHybrid()
    rates = ReactionRatesHybridJump()
    catalyst = CatalystStateDynamicJump()
    jacobian = JacobianLumpedSeidelJump()
    prefaktor = []
end
@kwdef struct BusscheJump <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterBussche()
    rates = ReactionRatesBusscheJump()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
@kwdef struct NestlerJump <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterNestler()
    rates = ReactionRatesNestlerJump()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
@kwdef struct BrilmanJump <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterBrilman()
    rates = ReactionRatesBrilmanJump()
    catalyst = CatalystStateNone()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
@kwdef struct AutoCatJump <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 1; -1 0 -1; -2 -3 1; 0 1 -1; 0 0 0]
    parameter = KineticParameterAutoCat()
    rates = ReactionRatesAutoCatJump()
    catalyst = CatalystStateDynamic()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
@kwdef struct LinearJump <: AbstractKineticSetup
    stoich = [1 1 0; 0 -1 -1; -1 0 1; -2 -3 -1; 0 1 1; 0 0 0]
    parameter = KineticParameterLinear()
    rates = ReactionRatesLinearJump()
    catalyst = CatalystStateDynamic()
    jacobian = JacobianLumpedNone()
    prefaktor = PrefaktorNone()
end
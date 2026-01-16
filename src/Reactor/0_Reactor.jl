"""
    struct Reactor(p)

Create a struct for reactor parameter and methods. p set the method modes:

    SetupVollbrecht()   
    SetupBerty()
    SetupBifurcationAnalysis()
    SetupBifurcationAnalysisMethane()
    SetupFpoCstr()
    SetupFpoPftr()
    SetupBertyNMPC()

field:

    properties  ::Any         #(reactorproperties defined by setup function)
    kinetic     ::Any         #(Kinetic defined by setup function)
    methods     ::any         #(defined by setup funktion)
"""
struct Reactor <: AbstractReactor
    properties
    kinetic
    methods
end
function Reactor(p::T) where {T<:AbstractExperimental}
    Reactor(
        ReactorProperties(p),
        Kinetic(p.kinetic),
        Experiments(p.xin_mess, p.xout_mess, nothing, p.files, p.idx, ObjectiveSS(), ObjectiveDyn(), p.subset),
    )
end
function Reactor(p::T) where {T<:AbstractBifurcation}
    Reactor(
        ReactorProperties(p.properties),
        Kinetic(p.kinetic),
        Bifurcation(p.predictor, p.corrector, nothing, nothing),
    )
end
function Reactor(p::T) where {T<:AbstractFPO}
    nall = p.n + (p.ordn - 2) * (p.n - 1)
    mall = p.m + (p.ordm - 2) * (p.m - 1)
    nc = 9
    Reactor(
        ReactorProperties(p.properties),
        Kinetic(p.kinetic),
        FPO(p.obj1, p.obj2, p.isothermal, p.n, p.ordn, p.m, p.ordm, SolutionFPO(nall, mall, p.numP, nc), nall, mall, p.setMfcCons, p.inert, p.setObjective, p.signal, p.bounds, p.predictor, p.optimize, p.N2, p.numP, p.optimize_temperature, p.inputs),
    )
end
function Reactor(p::T) where {T<:AbstractMPC}
    nall = p.n + (p.ordn - 2) * (p.n - 1)
    prop = ReactorProperties(p.properties)
    splco2 = Spline1D(prop.yin[2] * prop.ndotin)
    splco = Spline1D(prop.yin[3] * prop.ndotin)
    Reactor(
        prop,
        Kinetic(p.kinetic),
        MPC(p.obj1, p.obj2, p.n, p.ordn, p.sol, nall, p.init, p.horizon, (co2=splco2, co=splco), p.disturb, p.para, p.setpoint, p.err_int, p.deadtime)
    )
end

# Setup Reactor properties  
@kwdef struct SetupVollbrecht <: AbstractExperimental
    xin_mess = nothing
    xout_mess = nothing
    kinetic = Seidel()
    files = nothing
    idx = 20
    subset = 1:140
end
@kwdef struct SetupBerty <: AbstractExperimental
    xin_mess = nothing
    xout_mess = nothing
    kinetic = Seidel()
    files = nothing
    idx = 20
    subset = 1:50
end
@kwdef struct SetupNestler <: AbstractExperimental
    xin_mess = nothing
    xout_mess = nothing
    kinetic = Seidel()
    files = nothing
    idx = 20
    subset = 1:50
end


###Setup Bifurcation
@kwdef struct SetupBifurcationAnalysis <: AbstractBifurcation
    kinetic = Seidel()
    properties = SetupVollbrecht()
    predictor = nothing
    corrector = nothing
end
@kwdef struct SetupBifurcationAnalysisMethane <: AbstractBifurcation
    kinetic = Koschany()
    properties = SetupMethanation()
    predictor = nothing
    corrector = nothing
end

###Setup FPO###
@kwdef struct SetupFpoCstr <: AbstractFPO
    kinetic = SeidelJump()
    properties = SetupBerty()
    isothermal = true
    obj1 = ObjectiveFlowJump()
    obj2 = ObjectiveYieldJump()
    n = 25
    ordn = 3
    m = 1
    ordm = 2
    sol = nothing
    setMfcCons = false
    inert = true
    setObjective = 1
    signal = x -> 0.0
    bounds = (umax=[1, 1, 1, 1.0, 1, 1, 1, 1.1, 1, 753.15, 753.15, 10.0, 1.0, 1.0, 1.0], umin=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.005, 0.0, 0.00, 300, 473.15, 0.0, 0.0, 0.0, 0.0])
    predictor = PredictorFPO()
    optimize = :CO_N2_NIN
    N2 = 0.15
    numP = 40
    optimize_temperature = false
    inputs = zeros(16)
end

@kwdef struct SetupFpoPftr <: AbstractFPO
    kinetic = SeidelJump()
    properties = SetupNestler()
    isothermal = true
    obj1 = ObjectiveFlowJump()
    obj2 = ObjectiveYieldJump()
    n = 25
    ordn = 3
    m = 25
    ordm = 3
    sol = nothing
    setMfcCons = false
    inert = true
    setObjective = 1
    signal = x -> 0.0
    bounds = (umax=[1, 1, 1, 1.0, 1, 1, 1.0, 1.1, 1, 753.15, 753.15, 10.0, 1.0, 1.0, 1.0], umin=[0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.05, -1.1, 0.00, 300, 453.15, 0.0, 0.0, 0.0, 0.0])
    predictor = PredictorFPO()
    optimize = :CO_N2_NIN
    N2 = 0.15
    numP = 40
    optimize_temperature = false
    inputs = zeros(16)
end

###Setup NMPC###

@kwdef struct SetupMPCJump <: AbstractMPC
    kinetic = SeidelJump()
    properties = SetupBerty()
    obj1 = nothing
    obj2 = nothing
    n = 100
    ordn = 3
    sol = nothing
    init = nothing
    horizon = 600
    input = nothing
    disturb = nothing
    para = ones(3)
    setpoint = 2
    err_int = 0
    deadtime = 120
end
@kwdef struct SetupMPCPlant <: AbstractMPC
    kinetic = Seidel()
    properties = SetupBerty()
    obj1 = nothing
    obj2 = nothing
    n = 100
    ordn = 3
    sol = nothing
    init = nothing
    horizon = 600
    input = nothing
    disturb = nothing
    para = ones(3)
    setpoint = 2
    err_int = 0
    deadtime = 120
end

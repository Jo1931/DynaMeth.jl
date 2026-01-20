using DynaMeth
using LaTeXStrings
using JuMP
using Plots

import HSL                  # Wrapper-Package
HSL.LIBHSL_isfunctional()   # → true, if installed

###########Steady State Front################
#Setup Reactor
fpo = SetupFpoPftr(n=2, ordn=2, m=2, ordm=2, numP=40, N2=0.15, kinetic=SeidelJump(), optimize_temperature=false, isothermal=true)
reacSS = Reactor(fpo)
reacSS.properties.Tin = 503.15
reacSS.properties.P = 60 * 1e5
reacSS.methods.signal = x -> 0

#Left Bound
reacSS.methods.setObjective = 1
init_reactor = reacSS.methods.predictor(reacSS, 1)
optimizeRPLUG!(init_reactor, 0.0, 1, reacSS, sca=1000.0, linear_solver="ma97");

#Right Bound
reacSS.methods.setObjective = 2
init_reactor = reacSS.methods.predictor(reacSS, 2)
optimizeRPLUG!(init_reactor, 0.0, reacSS.methods.numP, reacSS; sca=1.0, linear_solver="ma97");

#ϵ-constraint
epsilon = LinRange(reacSS.methods.sol.obj2[1], reacSS.methods.sol.obj2[end], reacSS.methods.numP)

## Front ##
reacSS.methods.setObjective = 1
for i = reacSS.methods.numP-1:-1:2
    init_reactor_s = reacSS.methods.predictor(reacSS, i, height=0.3, direction=:forward)
    optimizeRPLUG!(init_reactor_s, epsilon[i], i, reacSS; sca=1.0, linear_solver="ma97")
    print(string(i) * " ")
end

###########Periodic Front####################
#Setup Reactor
fpo = SetupFpoPftr(n=50, ordn=2, m=2, ordm=2, numP=40, N2=0.15, kinetic=SeidelJump(), optimize_temperature=false, isothermal=true)
reac = Reactor(fpo)
reac.properties.Tin = 503.15
reac.properties.P = 60 * 1e5
reac.methods.signal = x -> sig(x)

#Left Bound
reac.methods.setObjective = 1
init_reactor = reac.methods.predictor(reac, 1)
optimizeRPLUG!(init_reactor, 0.0, 1, reac, sca=10000.0, linear_solver="ma97");

#Right Bound
reac.methods.setObjective = 2
init_reactor = reacSS.methods.predictor(reac, 1)
optimizeRPLUG!(init_reactor, 0.0, reac.methods.numP, reac; sca=0.01, linear_solver="ma97");

#ϵ-constraint
epsilon = LinRange(reac.methods.sol.obj2[1], reac.methods.sol.obj2[end], reac.methods.numP)

## Front ##
reac.methods.setObjective = 1
for i = reac.methods.numP-1:-1:2
    init_reactor_s = reac.methods.predictor(reac, i, height=0.1, direction=:forward)
    optimizeRPLUG!(init_reactor_s, epsilon[i], i, reac; sca=10000000.0, linear_solver="ma97")
    print(string(i) * " ")
end

###Plot###
plot(reacSS.methods.sol, c=:black, label="Steady State", framestyle=:box, ylabel=L"\Phi_1\mathrm{\ in \ mmol/min/kg}")
plot!(reac.methods.sol, c=:red, label="Forced Periodic Operation", xlabel=L"\Phi_2")


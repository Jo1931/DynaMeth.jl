using DynaMeth
using OptimizationOptimJL
using Plots

#read experimental data
xin, xout = createSteadyStates("Volllbrecht")
dynamic = DynamicVollbrecht()

#initialize Reactor (Nestler kinetics)
vollbrecht = SetupVollbrecht(xin_mess=xin, xout_mess=xout, kinetic=Nestler(), files=dynamic)
cstr = Reactor(vollbrecht)

#only experiments 1:79 (pure CO experiments excluded)
cstr.methods.subset = 1:79
p = plotReactorSteady(cstr)
cstr.methods.objectiveSS(cstr)

#run optimization
n = 10          # can be increased for higher accuracy
for i = 1:n
    solve_opt(cstr; iter=500)   #Nelder-Mead
    try
        solve_opt(cstr; solver=Optim.ParticleSwarm(; n_particles=5), iter=200)  #ParticleSwarm
    catch
    end
end

#Plot
print("Residual:" * string(cstr.methods.objectiveSS(cstr)))
p = plotReactorSteady!(p, cstr)

















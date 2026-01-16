struct ObjectiveDyn end
# function (res::ObjectiveDyn)(u, cstr::Reactor, ::SetupBerty)
#     files = cstr.experiments.files
#     setDecisionVar!(cstr, u)
#     t, x_out, sol, xin, err = simulateAllExperiments(cstr, files)
#     penalty = 0
#     for i = 7:14
#         penalty = penalty + minimum([0, cstr.kinetic.parameter.beta[i]]) .^ 2
#     end
#     for i = 1:2
#         penalty = penalty + minimum([0, cstr.kinetic.parameter.k[i]]) .^ 2
#     end
#     res = sum(err) / length(err) + penalty * 1e8
# end
# function (res::ObjectiveDyn)(u, cstr::Reactor)
#     setDecisionVar!(cstr, u)
#     data1 = cstr.experiments.files
#     ts = data1[!, "dynamic_1"] * 60
#     ts, x_out, sol, xin = simulateVollbrechtDyn(cstr, ts)
#     ch3oh = sol[1, :] - data1[!, "dynamic_2"]
#     co2 = sol[2, :] - data1[!, "dynamic_3"]
#     co = sol[3, :] - data1[!, "dynamic_4"]
#     out = sum(ch3oh .^ 2) + sum(co2 .^ 3) + sum(co .^ 4)
#     return out
# end
function (res::ObjectiveDyn)(x, cstr::T) where {T<:Union{Reactor,ComponentVector}}
    #x = unwrappedParam(u, cstr)
    #setDecisionVar!(cstr, x)

    u = unwrappedParam(x, cstr)

    #penalty = penaltyFunc(cstr)
    data1 = cstr.methods.files

    ts = data1[!, "dynamic_1"] * 60
    set_kinetic_param!(cstr, u)
    ts, x_out, sol, xin = simulateVollbrechtDyn(cstr, ts)

    ch3oh = ((sol[1, :] - data1[!, "dynamic_2"]) ./ (data1[!, "dynamic_2"] .+ 1e-1)) .^ 2
    co2 = ((sol[2, :] - data1[!, "dynamic_3"]) ./ (data1[!, "dynamic_3"] .+ 1e-1)) .^ 2
    co = ((sol[3, :] - data1[!, "dynamic_4"]) ./ (data1[!, "dynamic_3"] .+ 1e-1)) .^ 2
    dyn = sum(ch3oh .^ 2) + sum(co2 .^ 3) + sum(co .^ 4)
    res = dyn
end
# function (res::ObjectiveDyn)(u, cstr::Reactor)
#     res(u, cstr, cstr.properties.type)
# end

struct ObjectiveYieldJump end
function (f::ObjectiveYieldJump)(model, x, n_in2, n_in3, tspan, n, mcat)


    #calculate mean values
    @expression(model, obj2, model[:obj1] * mcat / 6e4 / (n_in2 + n_in3))
    nothing
end
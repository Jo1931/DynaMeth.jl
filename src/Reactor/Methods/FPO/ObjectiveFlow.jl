struct ObjectiveFlowJump end
function (res::ObjectiveFlowJump)(model, x, tspan, n, mcat)
    # Objective functions
    @expression(model, obj1, sum((x[9, j] * x[1, j] * 0.5 + x[9, j-1] * x[1, j-1] * 0.5) * (tspan[j] - tspan[j-1]) for j = 2:n) / mcat * 6e4)
    nothing
end



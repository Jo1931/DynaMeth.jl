function calcSysMatrix(cstr::Reactor, xn)

    f(x) = vcat(inv(calcLHSbif(vcat(x, xn[9]), cstr)[1:8, 1:8]) * calcRHSsingle(vcat(x, xn[9]), cstr)[1:8])
    ForwardDiff.jacobian(f, xn[1:8])
end


# function calcSysMatrixCascade(cstr::Reactor, xn)
#     nc = cstr.properties.num_casc
#     v = collect(1:9*nc)
#     deleteat!(v, 9:9:length(v))

#     function f(x)
#         xx = vcat(x[1:8], xn[9])
#         for i = 2:nc
#             h = vcat(x[8*(i-1)+1:8*i], xn[9*i])
#             xx = vcat(xx, h)
#         end
#         RHS = calcRHScascade(xx, cstr::Reactor)
#         RHSn = zeros(eltype(x), 8 * nc)
#         for i = 1:nc
#             LHS = calcLHSbif(xx[(9*i-8):(9*i)], cstr)
#             LHS[9, 9] = 2
#             RHSn[(8*i-7):(8*i)] = inv(LHS[1:9, 1:9])[1:8, 1:9] * RHS[(9*i-8):(9*i)]
#         end
#         return RHSn
#     end
#     #f(x) = vcat(inv(calcLHSbif(x, cstr)[1:8, 1:8]) * calcRHSsingle(x, cstr)[1:8], 0)

#     ForwardDiff.jacobian(f, xn[v])
#     #FiniteDiff.finite_difference_jacobian(f, xn[v])
# end

function calcSysMatrixCascade(cstr::Reactor, xn)
    nc = cstr.properties.num_casc

    function f(xx)
        x = xx
        RHS = calcRHScascade(xx, cstr::Reactor)
        RHSn = zeros(eltype(x), 9 * nc)
        for i = 1:nc
            LHS = calcLHSbif(xx[(9*i-8):(9*i)], cstr)
            LHS[9, 9] = 1
            RHSn[(9*i-8):(9*i)] = inv(LHS[1:9, 1:9]) * RHS[(9*i-8):(9*i)]
        end
        return RHSn
    end
    #f(x) = vcat(inv(calcLHSbif(x, cstr)[1:8, 1:8]) * calcRHSsingle(x, cstr)[1:8], 0)
    #print(sum(f(xn)))
    #print("  ")
    ForwardDiff.jacobian(f, xn)
    #FiniteDiff.finite_difference_jacobian(f, xn)
end

